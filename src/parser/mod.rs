/// Module to parse `hp.obo` file
pub mod hp_obo {

//! This will work the following way:
//! 1. Add new terms
//! 2. Replace obsolete ones
//! 3. Connect parents

}

/// Module to parse HPO - Gene associations from `phenotype_to_genes.txt` file
pub mod phenotype_to_genes {
    use std::fs::File;
    use std::io::BufRead;
    use std::io::BufReader;

    use crate::HpoTermId;
    use crate::Ontology;

    /// Quick and dirty parser for development and debugging
    pub fn parse(file: &str, collection: &mut Ontology) {
        let file = File::open(file).unwrap();
        let reader = BufReader::new(file);
        for line in reader.lines() {
            let line = line.unwrap();
            if line.starts_with('#') {
                continue;
            }
            let cols: Vec<&str> = line.trim().split('\t').collect();
            let gene_id = collection.add_gene(cols[3], cols[2]).unwrap();
            let term_id = HpoTermId::from(cols[0]);
            collection.link_gene_term(&term_id, gene_id);

            collection
                .gene_mut(&gene_id)
                .expect("Cannot find gene")
                .add_term(term_id);
        }
    }
}

/// Module to parse HPO - OmimDisease associations from `phenotype.hpoa` file
pub mod phenotype_hpoa {
    use crate::HpoTermId;
    use std::fs::File;
    use std::io::BufRead;
    use std::io::BufReader;

    use crate::Ontology;

    struct Omim<'a> {
        id: &'a str,
        name: &'a str,
        hpo_id: HpoTermId,
    }

    fn parse_line(line: &str) -> Option<Omim<'_>> {
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

        let (_, omim_id) = cols[0].split_once(':').unwrap();

        Some(Omim {
            id: omim_id,
            name: cols[1],
            hpo_id: HpoTermId::from(cols[3]),
        })
    }

    /// Quick and dirty parser for development and debugging
    pub fn parse(file: &str, collection: &mut Ontology) {
        let file = File::open(file).unwrap();
        let reader = BufReader::new(file);
        for line in reader.lines() {
            let line = line.unwrap();
            if let Some(omim) = parse_line(&line) {
                let omim_disease_id = collection.add_omim_disease(omim.name, omim.id).unwrap();
                collection.link_omim_disease_term(&omim.hpo_id, omim_disease_id);

                collection
                    .omim_disease_mut(&omim_disease_id)
                    .expect("Cannot find gene")
                    .add_term(omim.hpo_id);
            }
        }
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
