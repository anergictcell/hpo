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

    let mut genes: Vec<&str> = Vec::new();
    for term in collection.hpos() {
        genes.clear();
        for gene in term.genes() {
            genes.push(gene.name());
        }
        genes.sort();
        println!("{}\t{}", term.id(), genes.join(","));
    }
}

/*
Python comparison code
from pyhpo import Ontology

_ = Ontology()

with open("term2gene.py.txt", "w") as fh:
    for term in Ontology:
        genes = []
        for gene in term.genes:
            genes.append(gene.name)
        genes = ",".join(sorted(genes))
        _ = fh.write(f"{term.id}\t{genes}\n")

*/
