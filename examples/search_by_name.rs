//! Prints every term and its associated genes

use hpo::Ontology;

fn main() {
    let ontology = Ontology::from_binary("tests/ontology.hpo").unwrap();
    let cystinosis = ontology.disease_by_name("Cystinosis").unwrap();
    println!("first match: {:?}", cystinosis.name());
    for result in ontology.diseases_by_name("Cystinosis") {
        println!("{:?}", result.name());
    }
}
