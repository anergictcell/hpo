//! Lists all common ancestors of two HPO terms, along with metadata
//! like number of associated OMIM diseases or genes and the OMIM-based
//! Information Content

use hpo::{HpoTermId, Ontology};

fn common_ancestors(termid1: HpoTermId, termid2: HpoTermId, ontology: &Ontology) {
    let term1 = ontology.hpo(termid1).unwrap();
    let term2 = ontology.hpo(termid2).unwrap();
    for term in &term1.common_ancestors(&term2) {
        println!(
            "Term {} | IC (Omim) {} | IC (Gene) {} | nOmim {} | nGene {}",
            term.id(),
            term.information_content().omim_disease(),
            term.information_content().gene(),
            term.omim_diseases().count(),
            term.genes().count()
        );
    }
}

fn main() {
    let ontology = Ontology::from_binary("tests/ontology.hpo").unwrap();

    let mut args = std::env::args();
    if args.len() == 3 {
        let term_id1 = args.nth(1).unwrap();
        let term_id2 = args.next().unwrap();
        common_ancestors(term_id1.into(), term_id2.into(), &ontology);
    } else {
        println!("Show all common ancestors of 2 HPO terms\n\n");
        println!("Usage\ncommon_ancestors <HPO-TERM> <HPO-TERM>");
        println!("e.g.:\ncommon_ancestors HP:0007768 HP:0000631\n");
    }
}
