//! Prints every term and its associated genes

use hpo::annotations::Disease;
use hpo::Ontology;

fn main() {
    let ontology = Ontology::from_binary("tests/ontology.hpo").unwrap();
    let congenital = ontology.omim_disease_by_name("congenital").unwrap();
    println!("first match: {:?}", congenital.name());
    for result in ontology.omim_diseases_by_name("congenital") {
        println!("{:?}", result.name());
    }
}
