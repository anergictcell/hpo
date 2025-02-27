//! Parsing the HPO master data provided by Jax

use std::path::Path;

use builder::Builder;

use crate::{ontology::builder, HpoResult, Ontology};

pub(crate) mod binary;
/// Module to parse `hp.obo` file
pub(crate) mod hp_obo;

/// Module to parse HPO - `Gene` associations
///
/// It contains functions to parse `genes_to_phenotype.txt` and
/// `phenotype_to_genes.txt` input files
pub(crate) mod gene_to_hpo {

    use crate::annotations::GeneId;
    use crate::ontology::builder::ConnectedTerms;
    use crate::ontology::Builder;
    use crate::parser::Path;
    use crate::HpoError;
    use crate::HpoResult;
    use std::fs::File;
    use std::io::BufRead;
    use std::io::BufReader;

    use crate::HpoTermId;

    struct ParsedGene<'a> {
        ncbi_id: GeneId,
        symbol: &'a str,
        hpo: HpoTermId,
    }

    impl<'a> ParsedGene<'a> {
        fn try_new(ncbi_id: &'a str, symbol: &'a str, hpo: &'a str) -> HpoResult<Self> {
            let hpo = HpoTermId::try_from(hpo)?;
            let ncbi_id = GeneId::try_from(ncbi_id)?;
            Ok(Self {
                ncbi_id,
                symbol,
                hpo,
            })
        }
    }

    /// Removes the first (header) line.
    ///
    /// TODO: Update this once <https://doc.rust-lang.org/std/io/trait.BufRead.html#method.skip_until> is stable
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
    ///
    /// ```text
    /// 10  NAT2    HP:0000007  Autosomal recessive inheritance         -       OMIM:243400
    /// ```
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

        // Column 2 is the HPO-name, which we don't need
        if cols.next().is_none() {
            return Err(HpoError::InvalidInput(line.to_string()));
        };

        // Column 3 is the NCBI-ID of the gene
        let Some(ncbi_id) = cols.next() else {
            return Err(HpoError::InvalidInput(line.to_string()));
        };

        // Column 4 is the gene symbol
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
        builder: &mut Builder<ConnectedTerms>,
    ) -> HpoResult<()> {
        parse(file, builder, genes_to_phenotype_line)
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
        builder: &mut Builder<ConnectedTerms>,
    ) -> HpoResult<()> {
        parse(file, builder, phenotype_to_gene_line)
    }

    /// Parses a file to connect genes to HPO terms
    fn parse<P: AsRef<Path>, F: Fn(&str) -> HpoResult<ParsedGene<'_>>>(
        file: P,
        builder: &mut Builder<ConnectedTerms>,
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
            builder.annotate_gene(gene.ncbi_id, gene.symbol, gene.hpo)?;
        }
        Ok(())
    }

    #[cfg(test)]
    mod test {
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
        fn test_remove_header_hpo_id() {
            let x = "hpo_id\txyz\n10\tNAT2\n".as_bytes();
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
    }

    #[cfg(test)]
    mod test_genes_to_phenotype {
        use crate::annotations::AnnotationId;

        use super::*;

        #[test]
        fn test_parse_correct_line_with_newline() {
            let line = "10\tNAT2\tHP:0000007\tfoobar\n";
            let res = genes_to_phenotype_line(line).expect("This line should parse correctly");
            assert_eq!(res.ncbi_id.as_u32(), 10);
            assert_eq!(res.symbol, "NAT2");
            assert_eq!(res.hpo.as_u32(), 7u32);
        }

        #[test]
        fn test_parse_correct_line_without_newline() {
            let line = "10\tNAT2\tHP:0000007\tfoobar";
            let res = genes_to_phenotype_line(line).expect("This line should parse correctly");
            assert_eq!(res.ncbi_id.as_u32(), 10);
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

    #[cfg(test)]
    mod test_phenotpye_to_genes {
        use super::*;
        use crate::annotations::AnnotationId;

        #[test]
        fn test_parse_correct_line_with_newline() {
            let line = "HP:0000007\tAbnormality of body height\t10\tNAT2\tfoobar\n";
            let res = phenotype_to_gene_line(line).expect("This line should parse correctly");
            assert_eq!(res.ncbi_id.as_u32(), 10);
            assert_eq!(res.symbol, "NAT2");
            assert_eq!(res.hpo.as_u32(), 7u32);
        }

        #[test]
        fn test_parse_correct_line_without_newline() {
            let line = "HP:0000007\tAbnormality of body height\t10\tNAT2\tfoobar";
            let res = phenotype_to_gene_line(line).expect("This line should parse correctly");
            assert_eq!(res.ncbi_id.as_u32(), 10);
            assert_eq!(res.symbol, "NAT2");
            assert_eq!(res.hpo.as_u32(), 7u32);
        }

        #[test]
        fn test_parse_missing_id() {
            let line = "HP:0000007\tAbnormality of body height\tNAT2\tfoobar\n";
            let res = phenotype_to_gene_line(line);
            assert!(res.is_err());
        }

        #[test]
        fn test_parse_empty_symbol() {
            let line = "HP:0000007\tAbnormality of body height\t10\t\tfoobar\n";
            let res = phenotype_to_gene_line(line).expect("This line should parse correctly");
            assert_eq!(res.ncbi_id.as_u32(), 10);
            assert_eq!(res.symbol, "");
            assert_eq!(res.hpo.as_u32(), 7u32);
        }

        #[test]
        fn test_parse_missing_hpo() {
            let line = "Abnormality of body height\t10\tNAT2\tfoobar\n";
            let res = phenotype_to_gene_line(line);
            assert!(res.is_err());
        }

        #[test]
        fn test_parse_invalid_hpo() {
            let line = "HP:0000007A\tAbnormality of body height\t10\tNAT2\tfoobar\n";
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
pub(crate) mod disease_to_hpo {
    use crate::annotations::OmimDiseaseId;
    use crate::annotations::OrphaDiseaseId;
    use crate::ontology::builder::ConnectedTerms;
    use crate::ontology::Builder;
    use crate::HpoError;
    use crate::HpoResult;
    use crate::HpoTermId;
    use std::fs::File;
    use std::io::BufRead;
    use std::io::BufReader;
    use std::path::Path;

    enum DiseaseKind<'a> {
        Omim(DiseaseComponents<'a>),
        Orpha(DiseaseComponents<'a>),
    }

    struct DiseaseComponents<'a> {
        id: &'a str,
        name: &'a str,
        hpo_id: HpoTermId,
    }

    impl DiseaseComponents<'_> {
        fn omim_disease_id(&self) -> HpoResult<OmimDiseaseId> {
            OmimDiseaseId::try_from(self.id)
        }

        fn orpha_disease_id(&self) -> HpoResult<OrphaDiseaseId> {
            OrphaDiseaseId::try_from(self.id)
        }
    }

    fn parse_line(line: &str) -> HpoResult<Option<DiseaseKind<'_>>> {
        if line.starts_with("OMIM") {
            Ok(parse_disease_components(line)?.map(DiseaseKind::Omim))
        } else if line.starts_with("ORPHA") {
            Ok(parse_disease_components(line)?.map(DiseaseKind::Orpha))
        } else {
            Ok(None)
        }
    }

    fn parse_disease_components(line: &str) -> HpoResult<Option<DiseaseComponents>> {
        let mut cols = line.trim().splitn(5, '\t');

        let Some(id_col) = cols.next() else {
            return Err(HpoError::InvalidInput(line.to_string()));
        };

        let Some((_, disease_id)) = id_col.split_once(':') else {
            return Err(HpoError::InvalidInput(line.to_string()));
        };

        let Some(disease_name) = cols.next() else {
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

        Ok(Some(DiseaseComponents {
            id: disease_id,
            name: disease_name,
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
    pub fn parse<P: AsRef<Path>>(file: P, builder: &mut Builder<ConnectedTerms>) -> HpoResult<()> {
        let filename = file.as_ref().display().to_string();
        let file = File::open(file).map_err(|_| HpoError::CannotOpenFile(filename))?;
        let reader = BufReader::new(file);
        for line in reader.lines() {
            let line = line.unwrap();
            match parse_line(&line)? {
                Some(DiseaseKind::Omim(omim)) => {
                    builder.annotate_omim_disease(
                        omim.omim_disease_id()?,
                        omim.name,
                        omim.hpo_id,
                    )?;
                }
                Some(DiseaseKind::Orpha(orpha)) => {
                    builder.annotate_orpha_disease(
                        orpha.orpha_disease_id()?,
                        orpha.name,
                        orpha.hpo_id,
                    )?;
                }
                _ => {}
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
        fn test_correct_orpha() {
            let s = "ORPHA:600171\tGonadal agenesis\t\tHP:0000055\tOMIM:600171\tTAS\tP\tHPO:skoehler[2014-11-27]";
            let orpha = parse_line(s)
                .expect("This line has the correct format")
                .expect("Line describes an Omim disease");
            if let DiseaseKind::Orpha(orpha) = orpha {
                assert_eq!(orpha.name, "Gonadal agenesis");
                assert_eq!(orpha.id, "600171");
                assert_eq!(orpha.hpo_id, "HP:0000055");
            } else {
                panic!("Orpha line should be parsed as Orpha correctly");
            }
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
            if let DiseaseKind::Omim(omim) = omim {
                assert_eq!(omim.name, "Gonadal agenesis");
                assert_eq!(omim.id, "600171");
                assert_eq!(omim.hpo_id, "HP:0000055");
            } else {
                panic!("Omim line should be parsed as Omim correctly");
            }
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
) -> HpoResult<Ontology> {
    let builder = Builder::new();
    let builder = hp_obo::read_obo_file(obo_file, builder)?;
    let mut builder = builder.connect_all_terms();
    gene_to_hpo::parse_phenotype_to_genes(gene_file, &mut builder)?;
    disease_to_hpo::parse(disease_file, &mut builder)?;
    builder
        .calculate_information_content()?
        .build_with_defaults()
}

pub(crate) fn load_from_jax_files<P: AsRef<Path>>(
    obo_file: P,
    gene_file: P,
    disease_file: P,
) -> HpoResult<Ontology> {
    let builder = Builder::new();
    let builder = hp_obo::read_obo_file(obo_file, builder)?;
    let mut builder = builder.connect_all_terms();
    gene_to_hpo::parse_genes_to_phenotype(gene_file, &mut builder)?;
    disease_to_hpo::parse(disease_file, &mut builder)?;
    builder
        .calculate_information_content()?
        .build_with_defaults()
}
