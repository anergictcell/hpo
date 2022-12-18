use std::{fs::File, io::Write};

use hpo::Ontology;

fn main() {
    let ontology = Ontology::from_standard("./example_data/");
    println!("Ontology with {} terms", ontology.len());
    let mut args = std::env::args();
    if args.len() == 2 {
        let filename = args.nth(1).unwrap();
        let mut fh = File::create(filename).unwrap();
        match fh.write_all(&ontology.as_bytes()) {
            Ok(_) => println!("Saved output"),
            Err(err) => println!("Error: {}", err),
        };
    } else {
        println!("Please specify an output file")
    }
}
