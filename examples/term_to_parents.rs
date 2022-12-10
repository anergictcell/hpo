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
        collection.add_parent(cols[1].try_into().unwrap(), cols[0].try_into().unwrap());
    }

    collection.create_cache();

    parser::phenotype_to_genes::parse("phenotype_to_genes.txt", collection);
}

fn main() {
    let mut collection = Ontology::default();
    from_file(&mut collection);

    let mut parents: Vec<String> = Vec::new();
    for term in collection.hpos() {
        parents.clear();
        for parent in term.all_parents() {
            parents.push(parent.id().to_string());
        }
        parents.sort();
        println!("{}\t{}", term.id(), parents.join(","));
    }
}

/*
Python comparison code
from pyhpo import Ontology

_ = Ontology()

with open("term2parents.py.txt", "w") as fh:
    for term in Ontology:
        parents = []
        for parent in term.all_parents:
            parents.append(parent.id)
        parents = ",".join(sorted(parents))
        _ = fh.write(f"{term.id}\t{parents}\n")

*/
