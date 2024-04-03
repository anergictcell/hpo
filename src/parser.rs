//! Parsing the HPO master data provided by Jax

use std::path::Path;

use crate::{HpoResult, Ontology};

pub(crate) mod binary;
/// Module to parse `hp.obo` file
pub(crate) mod hp_obo;

/// Module to parse HPO - `Gene` associations
///
pub(crate) mod gene_to_hpo {
    use crate::parser::Path;
    use crate::HpoError;
    use crate::HpoResult;
    use std::fs::File;
    use std::io::BufRead;
    use std::io::BufReader;

    use crate::HpoTermId;
    use crate::Ontology;

    struct ParsedGene<'a> {
        ncbi_id: &'a str,
        symbol: &'a str,
        hpo: HpoTermId,
    }

    impl<'a> ParsedGene<'a> {
        fn try_new(ncbi_id: &'a str, symbol: &'a str, hpo: &'a str) -> HpoResult<Self> {
            let hpo = HpoTermId::try_from(hpo)?;
            Ok(Self {
                ncbi_id,
                symbol,
                hpo,
            })
        }
    }

    // Removes the first (header) line.
    // TODO: Update this once https://doc.rust-lang.org/std/io/trait.BufRead.html#method.skip_until is stable
    fn remove_header<R: BufRead>(reader: &mut R) -> HpoResult<()> {
        let mut trash = String::with_capacity(80);
        reader.read_line(&mut trash).map_err(|_| {
            HpoError::InvalidInput("Invalid data in genes_to_phenotype.txt".to_string())
        })?;
        if !trash.starts_with('#')
            && !trash.starts_with("ncbi_gene_id")
            && !trash.starts_with("hpo_id")
        {
            return Err(HpoError::InvalidInput(
                "genes_to_phenotype.txt file must contain a header".to_string(),
            ));
        }
        Ok(())
    }

    /// Parses a single line of `genes_to_phenotype.txt`
    ///
    /// and returns a `ParsedGene` struct with gene and HPO info
    fn genes_to_phenotype_line(line: &str) -> HpoResult<ParsedGene<'_>> {
        let mut cols = line.split('\t');

        // Column 1 is the NCBI-ID of the gene
        let Some(ncbi_id) = cols.next() else {
            return Err(HpoError::InvalidInput(line.to_string()));
        };

        // Column 2 is the gene symbol
        let Some(symbol) = cols.next() else {
            return Err(HpoError::InvalidInput(line.to_string()));
        };

        // Column 3 is the Hpo Term ID
        let Some(hpo) = cols.next() else {
            return Err(HpoError::InvalidInput(line.to_string()));
        };

        ParsedGene::try_new(ncbi_id, symbol, hpo)
    }

    /// Parses a single line of `phenotype_to_genes.txt`
    ///
    /// ```text
    /// HP:0000002  Abnormality of body height  81848   SPRY4       orphadata   ORPHA:432
    /// ```
    fn phenotype_to_gene_line(line: &str) -> HpoResult<ParsedGene<'_>> {
        let mut cols = line.split('\t');

        // Column 1 is the Hpo Term ID
        let Some(hpo) = cols.next() else {
            return Err(HpoError::InvalidInput(line.to_string()));
        };

        if cols.next().is_none() {
            return Err(HpoError::InvalidInput(line.to_string()));
        };

        // Column 1 is the NCBI-ID of the gene
        let Some(ncbi_id) = cols.next() else {
            return Err(HpoError::InvalidInput(line.to_string()));
        };

        // Column 2 is the gene symbol
        let Some(symbol) = cols.next() else {
            return Err(HpoError::InvalidInput(line.to_string()));
        };

        ParsedGene::try_new(ncbi_id, symbol, hpo)
    }

    /// Parse `genes_to_phenotype.txt` file
    ///
    /// ```text
    /// ncbi_gene_id    gene_symbol hpo_id  hpo_name    frequency   disease_id
    /// 10  NAT2    HP:0000007  Autosomal recessive inheritance         -       OMIM:243400
    /// 10  NAT2    HP:0001939  Abnormality of metabolism/homeostasis   -       OMIM:243400
    /// 16  AARS1   HP:0002460  Distal muscle weakness                  15/15   OMIM:613287
    /// ```
    pub fn parse_genes_to_phenotype<P: AsRef<Path>>(
        file: P,
        ontology: &mut Ontology,
    ) -> HpoResult<()> {
        parse(file, ontology, genes_to_phenotype_line)
    }

    /// Parse `phenotype_to_genes.txt` file
    ///
    /// ```text
    /// #Format: HPO-id<tab>HPO label<tab>entrez-gene-id<tab>entrez-gene-symbol<tab>Additional Info from G-D source<tab>G-D source<tab>disease-ID for link
    /// HP:0000002  Abnormality of body height  81848   SPRY4       orphadata   ORPHA:432
    /// HP:0000002  Abnormality of body height  204219  CERS3       orphadata   ORPHA:79394
    /// HP:0000002  Abnormality of body height  51360   MBTPS2  -   mim2gene    OMIM:308205
    /// ```
    pub fn parse_phenotype_to_genes<P: AsRef<Path>>(
        file: P,
        ontology: &mut Ontology,
    ) -> HpoResult<()> {
        parse(file, ontology, phenotype_to_gene_line)
    }

    /// Parses a file to connect genes to HPO terms
    fn parse<P: AsRef<Path>, F: Fn(&str) -> HpoResult<ParsedGene<'_>>>(
        file: P,
        ontology: &mut Ontology,
        parse_line: F,
    ) -> HpoResult<()> {
        let filename = file.as_ref().display().to_string();
        let file = File::open(file).map_err(|_| HpoError::CannotOpenFile(filename))?;
        let mut reader = BufReader::new(file);

        remove_header(&mut reader)?;

        for line in reader.lines() {
            let line = line.map_err(|_| {
                HpoError::InvalidInput("Invalid data in genes_to_phenotype.txt".to_string())
            })?;

            let gene = parse_line(&line)?;

            let gene_id = ontology.add_gene(gene.symbol, gene.ncbi_id)?;
            ontology.link_gene_term(gene.hpo, gene_id)?;
            ontology
                .gene_mut(&gene_id)
                .expect("Gene is present because it was just add_omim_disease")
                .add_term(gene.hpo);
        }
        Ok(())
    }

    #[cfg(test)]
    mod test_gene_pheno_parsing {
        use crate::annotations::AnnotationId;

        use super::*;

        #[test]
        fn test_remove_header_ncbi_gene() {
            let x = "ncbi_gene_id\txyz\n10\tNAT2\n".as_bytes();
            let mut reader = BufReader::new(x);
            assert!(remove_header(&mut reader).is_ok());

            let mut lines = reader.lines();
            assert_eq!(lines.next().unwrap().unwrap(), "10\tNAT2");
            assert!(lines.next().is_none());
        }

        #[test]
        fn test_remove_header_hashtag() {
            let x = "#foobar\txyz\n10\tNAT2\n".as_bytes();
            let mut reader = BufReader::new(x);
            assert!(remove_header(&mut reader).is_ok());

            let mut lines = reader.lines();
            assert_eq!(lines.next().unwrap().unwrap(), "10\tNAT2");
            assert!(lines.next().is_none());
        }

        #[test]
        fn test_remove_header_fails() {
            let x = "foobar\txyz\n10\tNAT2\n".as_bytes();
            let mut reader = BufReader::new(x);
            assert!(remove_header(&mut reader).is_err());
        }

        #[test]
        fn test_parse_correct_line() {
            let line = "10\tNAT2\tHP:0000007\tfoobar\n";
            let res = genes_to_phenotype_line(line).expect("This line should parse correctly");
            assert_eq!(res.ncbi_id, "10");
            assert_eq!(res.symbol, "NAT2");
            assert_eq!(res.hpo.as_u32(), 7u32);
        }

        #[test]
        fn test_parse_correct_line2() {
            let line = "10\tNAT2\tHP:0000007\tfoobar";
            let res = genes_to_phenotype_line(line).expect("This line should parse correctly");
            assert_eq!(res.ncbi_id, "10");
            assert_eq!(res.symbol, "NAT2");
            assert_eq!(res.hpo.as_u32(), 7u32);
        }

        #[test]
        fn test_parse_missing_id() {
            let line = "NAT2\tHP:0000007\tfoobar\n";
            let res = genes_to_phenotype_line(line);
            assert!(res.is_err());
        }

        #[test]
        fn test_parse_missing_symbol() {
            let line = "10\tHP:0000007\tfoobar\n";
            let res = genes_to_phenotype_line(line);
            assert!(res.is_err());
        }

        #[test]
        fn test_parse_missing_hpo() {
            let line = "10\tNAT2\tfoobar\n";
            let res = genes_to_phenotype_line(line);
            assert!(res.is_err());
        }

        #[test]
        fn test_parse_invalid_hpo() {
            let line = "10\tNAT2\tHP:000000A\tfoobar\n";
            let res = genes_to_phenotype_line(line);
            assert!(res.is_err());
        }
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
    gene_to_hpo::parse_phenotype_to_genes(gene_file, ontology)?;
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
    gene_to_hpo::parse_genes_to_phenotype(gene_file, ontology)?;
    phenotype_hpoa::parse(disease_file, ontology)?;
    Ok(())
}
