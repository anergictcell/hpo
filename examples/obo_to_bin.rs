use std::{fs::File, io::Write};

use hpo::Ontology;

fn main() {
    let mut args = std::env::args();
    if args.len() == 3 {
        let folder = args.nth(1).unwrap();
        let ontology = Ontology::from_standard(&folder).unwrap();
        println!("Ontology with {} terms", ontology.len());
        let filename = args.next().unwrap();
        let mut fh = File::create(filename).unwrap();
        match fh.write_all(&ontology.as_bytes()) {
            Ok(_) => println!("Saved output"),
            Err(err) => println!("Error: {}", err),
        };
    } else {
        println!("Please specify an output file")
    }
}
