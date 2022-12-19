use hpo::Ontology;

fn main() {
    let mut args = std::env::args();
    if args.len() == 2 {
        let filename = args.nth(1).unwrap();
        let ontology = Ontology::from_binary(filename).unwrap();
        println!("Ontology with {} terms", ontology.len());
    } else {
        println!("Please specify an output file")
    }
}
