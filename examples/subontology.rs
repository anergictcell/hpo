

use hpo::{Ontology};



fn main() {
    let ontology = Ontology::from_binary("tests/ontology.hpo").unwrap();
    let sub = ontology.sub_ontology(
        ontology.hpo(1u32.into()).expect("root doesn't exist"),
        vec![
            ontology.hpo(11017u32.into()).expect("11017 doesnt exist"),
            ontology.hpo(25454u32.into()).expect("25454"),
            ontology.hpo(12285u32.into()).expect("12285"),
            ontology.hpo(864u32.into()).expect("864"),
            ontology.hpo(10662u32.into()).expect("10662"),
            ontology.hpo(12638u32.into()).expect("12638"),
            
        ]
    ).unwrap();

    println!("{}", &sub.as_mermaid());
}
