//! Parsing the HPO master data provided by Jax

use std::path::Path;

use crate::{HpoResult, Ontology};

pub(crate) mod binary;
/// Module to parse `hp.obo` file
pub(crate) mod hp_obo;

/// Module to parse HPO - `Gene` associations from `genes_to_phenotype.txt` file
pub(crate) mod genes_to_phenotype {
    use crate::parser::Path;
    use crate::HpoError;
    use crate::HpoResult;
    use std::fs::File;
    use std::io::BufRead;
    use std::io::BufReader;

    use crate::HpoTermId;
    use crate::Ontology;

    /// Quick and dirty parser for development and debugging
    pub fn parse<P: AsRef<Path>>(file: P, ontology: &mut Ontology) -> HpoResult<()> {
        let filename = file.as_ref().display().to_string();
        let file = File::open(file).map_err(|_| HpoError::CannotOpenFile(filename))?;
        let reader = BufReader::new(file);
        for line in reader.lines() {
            let line = line.unwrap();
            // TODO: Check for the header outside of the `lines` iterator
            if line.starts_with('#') || line.starts_with("ncbi_gene_id") {
                continue;
            }
            let cols: Vec<&str> = line.trim().split('\t').collect();
            let gene_id = ontology.add_gene(cols[1], cols[0])?;
            let term_id = HpoTermId::try_from(cols[2])?;
            ontology.link_gene_term(term_id, gene_id)?;

            ontology
                .gene_mut(&gene_id)
                .expect("Cannot find gene {gene_id}")
                .add_term(term_id);
        }
        Ok(())
    }
}

/// Module to parse HPO - `Gene` associations from `phenotype_to_genes.txt` file
pub(crate) mod phenotype_to_genes {
    use crate::parser::Path;
    use crate::HpoError;
    use crate::HpoResult;
    use std::fs::File;
    use std::io::BufRead;
    use std::io::BufReader;

    use crate::HpoTermId;
    use crate::Ontology;

    /// Quick and dirty parser for development and debugging
    pub fn parse<P: AsRef<Path>>(file: P, ontology: &mut Ontology) -> HpoResult<()> {
        let filename = file.as_ref().display().to_string();
        let file = File::open(file).map_err(|_| HpoError::CannotOpenFile(filename))?;
        let reader = BufReader::new(file);
        for line in reader.lines() {
            let line = line.unwrap();
            // TODO: Check for the header outside of the `lines` iterator
            if line.starts_with('#') || line.starts_with("hpo_id") {
                continue;
            }
            let cols: Vec<&str> = line.trim().split('\t').collect();
            let gene_id = ontology.add_gene(cols[3], cols[2])?;
            let term_id = HpoTermId::try_from(cols[0])?;
            ontology.link_gene_term(term_id, gene_id)?;

            ontology
                .gene_mut(&gene_id)
                .expect("Cannot find gene {gene_id}")
                .add_term(term_id);
        }
        Ok(())
    }
}

/// Module to parse HPO - `OmimDisease` associations from `phenotype.hpoa` file
///
/// # Example line
///
/// ```text
/// OMIM:619340  Developmental and epileptic encephalopathy 96      HP:0011097  PMID:31675180  PCS  1/2  P  HPO:probinson[2021-06-21]
/// OMIM:609153  Pseudohyperkalemia                             NOT HP:0001878  PMID:2766660   PCS       P  HPO:lccarmody[2018-10-03]
/// ```
///
pub(crate) mod phenotype_hpoa {
    use crate::HpoError;
    use crate::HpoResult;
    use crate::HpoTermId;
    use std::fs::File;
    use std::io::BufRead;
    use std::io::BufReader;
    use std::path::Path;

    use crate::Ontology;

    struct Omim<'a> {
        id: &'a str,
        name: &'a str,
        hpo_id: HpoTermId,
    }

    fn parse_line(line: &str) -> HpoResult<Option<Omim<'_>>> {
        // TODO (nice to have): Add check to skip `database_id` header row
        // It is not strictly needed, because we're discarding non-OMIM rows
        if line.starts_with('#') {
            return Ok(None);
        }
        if !line.starts_with("OMIM") {
            return Ok(None);
        }

        let mut cols = line.trim().splitn(5, '\t');

        let Some(id_col) = cols.next() else {
            return Err(HpoError::InvalidInput(line.to_string()));
        };
        let Some((_, omim_id)) = id_col.split_once(':') else {
            return Err(HpoError::InvalidInput(line.to_string()));
        };

        let Some(omim_name) = cols.next() else {
            return Err(HpoError::InvalidInput(line.to_string()));
        };

        if let Some("NOT") = cols.next() {
            return Ok(None);
        };

        let hpo_id = if let Some(id) = cols.next() {
            HpoTermId::try_from(id)?
        } else {
            return Err(HpoError::InvalidInput(line.to_string()));
        };

        Ok(Some(Omim {
            id: omim_id,
            name: omim_name,
            hpo_id,
        }))
    }

    /// Quick and dirty parser for development and debugging
    ///
    /// # Errors
    ///
    /// - [`HpoError::CannotOpenFile`]: Source file not present or can't be opened
    /// - [`HpoError::ParseIntError`]: A line contains an invalid `omim_disease_id`
    /// - [`HpoError::DoesNotExist`]: A line contains a non-existing [`HpoTermId`]
    pub fn parse<P: AsRef<Path>>(file: P, ontology: &mut Ontology) -> HpoResult<()> {
        let filename = file.as_ref().display().to_string();
        let file = File::open(file).map_err(|_| HpoError::CannotOpenFile(filename))?;
        let reader = BufReader::new(file);
        for line in reader.lines() {
            let line = line.unwrap();
            if let Some(omim) = parse_line(&line)? {
                let omim_disease_id = ontology.add_omim_disease(omim.name, omim.id)?;
                ontology.link_omim_disease_term(omim.hpo_id, omim_disease_id)?;

                ontology
                    .omim_disease_mut(&omim_disease_id)
                    .expect("Omim disease was just added and cannot be missing")
                    .add_term(omim.hpo_id);
            }
        }
        Ok(())
    }

    #[cfg(test)]
    mod test_omim_parsing {
        use super::*;

        #[test]
        fn test_skip_comment() {
            let s = "#OMIM:600171\tGonadal agenesis\t\tHP:0000055\tOMIM:600171\tTAS\tP\tHPO:skoehler[2014-11-27]";
            assert!(parse_line(s)
                .expect("This line has the correct format")
                .is_none());
        }

        #[test]
        fn test_skip_not() {
            let s = "OMIM:600171\tGonadal agenesis\tNOT\tHP:0000055\tOMIM:600171\tTAS\tP\tHPO:skoehler[2014-11-27]";
            assert!(parse_line(s)
                .expect("This line has the correct format")
                .is_none());
        }

        #[test]
        fn test_skip_orpha() {
            let s = "ORPHA:600171\tGonadal agenesis\t\tHP:0000055\tOMIM:600171\tTAS\tP\tHPO:skoehler[2014-11-27]";
            assert!(parse_line(s)
                .expect("This line has the correct format")
                .is_none());
        }

        #[test]
        fn test_skip_orpha_not() {
            let s = "ORPHA:600171\tGonadal agenesis\tNOT\tHP:0000055\tOMIM:600171\tTAS\tP\tHPO:skoehler[2014-11-27]";
            assert!(parse_line(s)
                .expect("This line has the correct format")
                .is_none());
        }

        #[test]
        fn test_correct_omim() {
            let s = "OMIM:600171\tGonadal agenesis\t\tHP:0000055\tOMIM:600171\tTAS\tP\tHPO:skoehler[2014-11-27]";
            let omim = parse_line(s)
                .expect("This line has the correct format")
                .expect("Line describes an Omim disease");
            assert_eq!(omim.name, "Gonadal agenesis");
            assert_eq!(omim.id, "600171");
            assert_eq!(omim.hpo_id, "HP:0000055");
        }

        #[test]
        fn test_invalid_omim_id() {
            let s = "OMIM_600171\tGonadal agenesis\t\tHP:0000055\tOMIM:600171\tTAS\tP\tHPO:skoehler[2014-11-27]";
            assert!(parse_line(s).is_err());
        }

        #[test]
        fn test_invalid_hpo_id() {
            let s = "OMIM:600171\tGonadal agenesis\t\tH55\tOMIM:600171\tTAS\tP\tHPO:skoehler[2014-11-27]";
            assert!(parse_line(s).is_err());
        }

        #[test]
        fn test_invalid_input() {
            let s = "OMIM:600171 Gonadal agenesis  HP:0000055 OMIM:600171 TAS P HPO:skoehler[2014-11-27]";
            assert!(parse_line(s).is_err());
        }
    }
}

pub(crate) fn load_from_jax_files_with_transivitve_genes<P: AsRef<Path>>(
    obo_file: P,
    gene_file: P,
    disease_file: P,
    ontology: &mut Ontology,
) -> HpoResult<()> {
    hp_obo::read_obo_file(obo_file, ontology)?;
    phenotype_to_genes::parse(gene_file, ontology)?;
    phenotype_hpoa::parse(disease_file, ontology)?;
    Ok(())
}

pub(crate) fn load_from_jax_files<P: AsRef<Path>>(
    obo_file: P,
    gene_file: P,
    disease_file: P,
    ontology: &mut Ontology,
) -> HpoResult<()> {
    hp_obo::read_obo_file(obo_file, ontology)?;
    genes_to_phenotype::parse(gene_file, ontology)?;
    phenotype_hpoa::parse(disease_file, ontology)?;
    Ok(())
}
