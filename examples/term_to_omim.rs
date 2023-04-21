//! Prints every term and its associated Omim Disease

use hpo::Ontology;

fn main() {
    let ontology = Ontology::from_binary("tests/ontology.hpo").unwrap();

    let mut omim_diseases: Vec<String> = Vec::new();
    for term in ontology.hpos() {
        omim_diseases.clear();
        for disease in term.omim_diseases() {
            omim_diseases.push(disease.id().to_string());
        }
        omim_diseases.sort();
        println!("{}\t{}", term.id(), omim_diseases.join(","));
    }
}

/*
Python comparison code
from pyhpo import Ontology

_ = Ontology()

with open("term2omim.py.txt", "w") as fh:
    for term in Ontology:
        omim_diseases = []
        for disease in term.omim_diseases:
            omim_diseases.append(f"OMIM:{disease.id}")
        omim_diseases = ",".join(sorted(omim_diseases))
        _ = fh.write(f"{term.id}\t{omim_diseases}\n")

*/
