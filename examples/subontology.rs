use hpo::Ontology;
use std::fs::File;
use std::io::Write;

fn main() {
    let ontology = Ontology::from_binary("tests/ontology.hpo").unwrap();
    let sub = ontology
        .sub_ontology(
            ontology.hpo(1u32.into()).expect("root doesn't exist"),
            vec![
                ontology.hpo(11017u32.into()).expect("11017 doesnt exist"),
                ontology.hpo(25454u32.into()).expect("25454"),
                ontology.hpo(12285u32.into()).expect("12285"),
                ontology.hpo(864u32.into()).expect("864"),
                ontology.hpo(10662u32.into()).expect("10662"),
                ontology.hpo(12638u32.into()).expect("12638"),
                ontology.hpo(3581u32.into()).expect("3581"),
                ontology.hpo(7u32.into()).expect("7"),
                ontology.hpo(12648u32.into()).expect("12648"),
            ],
        )
        .unwrap();

    let mut args = std::env::args();
    if args.len() == 2 {
        let filename = args.nth(1).unwrap();
        let mut fh = File::create(filename).unwrap();
        match fh.write_all(&sub.as_bytes()) {
            Ok(_) => println!("Saved output"),
            Err(err) => println!("Error: {}", err),
        };
    }

    println!("{}", &sub.as_mermaid());
}
