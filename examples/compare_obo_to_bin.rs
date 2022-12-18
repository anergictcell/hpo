use std::time::SystemTime;

use hpo::annotations::{Gene, OmimDisease};
use hpo::Ontology;

fn main() {
    let mut args = std::env::args();
    if args.len() == 2 {
        let filename = args.nth(1).unwrap();
        let start = SystemTime::now();
        let bin = Ontology::from_binary(filename).unwrap();
        let end = SystemTime::now();
        let duration1 = end.duration_since(start).unwrap();

        let start = SystemTime::now();
        let obo = Ontology::from_standard("./example_data/");
        let end = SystemTime::now();
        let duration2 = end.duration_since(start).unwrap();

        println!(
            "Bin Ontology with {} terms in {} milliseconds",
            bin.len(),
            duration1.as_millis()
        );
        println!(
            "Obo Ontology with {} terms in {} milliseconds",
            obo.len(),
            duration2.as_millis()
        );

        for term_bin in bin.hpos() {
            let term_obo = obo
                .hpo(term_bin.id())
                .unwrap_or_else(|| panic!("Term {} not found", term_bin.id()));
            assert_eq!(term_bin.parents().count(), term_obo.parents().count());
            let s1 = term_bin.parent_ids();
            let s2 = term_obo.parent_ids();
            assert_eq!(s1.len(), s2.len());
            assert_eq!(s1.len(), (s1 & s2).len());

            let s1all = term_bin.all_parent_ids();
            let s2all = term_obo.all_parent_ids();
            assert_eq!(s1all.len(), s2all.len());
            assert_eq!(s1all.len(), (s1all & s2all).len());

            let s1all = term_bin.genes();
            let s2all = term_obo.genes();
            let mut t1genes: Vec<&Gene> = s1all.collect();
            let mut t2genes: Vec<&Gene> = s2all.collect();
            t1genes.sort_by(|ga, gb| ga.name().cmp(gb.name()));
            t2genes.sort_by(|ga, gb| ga.name().cmp(gb.name()));
            assert_eq!(t1genes, t2genes);

            let s1all = term_bin.omim_diseases();
            let s2all = term_obo.omim_diseases();
            let mut t1diseases: Vec<&OmimDisease> = s1all.collect();
            let mut t2diseases: Vec<&OmimDisease> = s2all.collect();
            t1diseases.sort_by(|ga, gb| ga.id().cmp(gb.id()));
            t2diseases.sort_by(|ga, gb| ga.id().cmp(gb.id()));
            assert_eq!(t1diseases, t2diseases);

            t1diseases.sort_by(|ga, gb| ga.name().cmp(gb.name()));
            t2diseases.sort_by(|ga, gb| ga.name().cmp(gb.name()));
            assert_eq!(t1diseases, t2diseases);

            assert_eq!(
                term_obo.information_content().gene(),
                term_bin.information_content().gene()
            );
            assert_eq!(
                term_obo.information_content().omim_disease(),
                term_bin.information_content().omim_disease()
            );
        }
        println!("Identical number of parents in each term");
        println!("Identical parents in each term");
        println!("Identical number of genes in each term");
        println!("Identical genes in each term");
        println!("Identical number of diseases in each term");
        println!("Identical diseases in each term");
        println!("Identical information content in each term");
    } else {
        println!("Please specify an output file")
    }
}
