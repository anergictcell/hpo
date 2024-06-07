//! Prints every term and its associated genes

use hpo::annotations::Disease;
use hpo::Ontology;

fn main() {
    let ontology = Ontology::from_binary("tests/ontology.hpo").unwrap();
    let cystinosis = ontology.omim_disease_by_name("Cystinosis").unwrap();
    println!("first match: {:?}", cystinosis.name());
    for result in ontology.omim_diseases_by_name("Cystinosis") {
        println!("{:?}", result.name());
    }
}
