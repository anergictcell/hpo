//! Parsing the HPO master data provided by Jax

use std::path::Path;

use crate::{HpoResult, Ontology};

pub(crate) mod binary;
/// Module to parse `hp.obo` file
pub(crate) mod hp_obo;

/// Module to parse HPO - `Gene` associations from `genes_to_phenotype.txt` file
pub(crate) mod genes_to_phenotype {
    use crate::parser::Path;
    use crate::HpoResult;
    use std::fs::File;
    use std::io::BufRead;
    use std::io::BufReader;

    use crate::HpoTermId;
    use crate::Ontology;

    /// Quick and dirty parser for development and debugging
    pub fn parse<P: AsRef<Path>>(file: P, ontology: &mut Ontology) -> HpoResult<()> {
        let file = File::open(file).unwrap();
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
    use crate::HpoResult;
    use std::fs::File;
    use std::io::BufRead;
    use std::io::BufReader;

    use crate::HpoTermId;
    use crate::Ontology;

    /// Quick and dirty parser for development and debugging
    pub fn parse<P: AsRef<Path>>(file: P, ontology: &mut Ontology) -> HpoResult<()> {
        let file = File::open(file).unwrap();
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
pub(crate) mod phenotype_hpoa {
    use crate::HpoError;
    use crate::HpoResult;
    use crate::HpoTermId;
    use std::fs::File;
    use std::io::BufRead;
    use std::io::BufReader;
    use std::path::Path;

    use tracing::error;

    use crate::Ontology;

    struct Omim<'a> {
        id: &'a str,
        name: &'a str,
        hpo_id: HpoTermId,
    }

    fn parse_line(line: &str) -> Option<Omim<'_>> {
        // TODO (nice to have): Add check to skip `database_id` header row
        // It is not strictly needed, because we're discarding non-OMIM rows
        if line.starts_with('#') {
            return None;
        }
        if !line.starts_with("OMIM") {
            return None;
        }

        let cols: Vec<&str> = line.trim().split('\t').collect();
        if cols[2] == "NOT" {
            return None;
        }

        let Some((_, omim_id)) = cols[0].split_once(':') else {
            error!("cannot parse OMIM ID from {}", cols[0]);
            return None;
        };

        let Ok(hpo_id) = HpoTermId::try_from(cols[3]) else {
            error!("invalid HPO ID: {}", cols[3]);
            return None;
        };

        Some(Omim {
            id: omim_id,
            name: cols[1],
            hpo_id,
        })
    }

    /// Quick and dirty parser for development and debugging
    pub fn parse<P: AsRef<Path>>(file: P, ontology: &mut Ontology) -> HpoResult<()> {
        let file = File::open(file).unwrap();
        let reader = BufReader::new(file);
        for line in reader.lines() {
            let line = line.unwrap();
            if let Some(omim) = parse_line(&line) {
                let omim_disease_id = ontology.add_omim_disease(omim.name, omim.id)?;
                ontology.link_omim_disease_term(omim.hpo_id, omim_disease_id)?;

                ontology
                    .omim_disease_mut(&omim_disease_id)
                    .ok_or(HpoError::DoesNotExist)?
                    .add_term(omim.hpo_id);
            }
        }
        Ok(())
    }

    #[cfg(test)]
    mod test_omim_parsing {
        use super::*;

        #[test]
        fn test_skip_not() {
            let s = "OMIM:600171\tGonadal agenesis\tNOT\tHP:0000055\tOMIM:600171\tTAS\tP\tHPO:skoehler[2014-11-27]";
            assert!(parse_line(s).is_none());
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
