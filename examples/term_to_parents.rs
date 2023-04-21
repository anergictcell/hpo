//! Prints every term and its parents

use hpo::Ontology;

fn main() {
    let ontology = Ontology::from_binary("tests/ontology.hpo").unwrap();

    let mut parents: Vec<String> = Vec::new();
    for term in ontology.hpos() {
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
