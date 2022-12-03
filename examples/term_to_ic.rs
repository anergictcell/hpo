use std::fs::File;
use std::io::BufRead;
use std::io::BufReader;

use hpo::parser;
use hpo::Ontology;

fn from_file(collection: &mut Ontology) {
    let file = File::open("terms.txt").unwrap();
    let reader = BufReader::new(file);
    for term in reader.lines() {
        collection.add_term_by_name(&term.unwrap());
    }

    let file = File::open("connections.txt").unwrap();
    let reader = BufReader::new(file);
    for line in reader.lines() {
        let line = line.unwrap();
        let cols: Vec<&str> = line.splitn(2, '\t').collect();
        collection.add_parent(cols[1].into(), cols[0].into());
    }

    collection.create_cache();

    parser::phenotype_to_genes::parse("phenotype_to_genes.txt", collection);
    parser::phenotype_hpoa::parse("phenotype.hpoa", collection);

    collection.calculate_information_content();
}

fn main() {
    let mut collection = Ontology::default();
    from_file(&mut collection);

    for term in collection.iter_terms() {
        println!(
            "{}\t{}\t{}",
            term.id(),
            term.information_content().gene(),
            term.information_content().omim_disease(),
        );
    }
}

/*
Python comparison code
from pyhpo import Ontology

_ = Ontology()

with open("term2ic.py.txt", "w") as fh:
    for term in Ontology:
        _ = fh.write(f"{term.id}\t{term.information_content.gene}\t{term.information_content.omim}\n")

*/
