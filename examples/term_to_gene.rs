//! Prints every term and its associated genes

use hpo::Ontology;

fn main() {
    let ontology = Ontology::from_binary("tests/ontology.hpo").unwrap();

    let mut genes: Vec<&str> = Vec::new();
    for term in ontology.hpos() {
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
