use std::{fs::File, io::Write};

use hpo::Ontology;

fn main() {
    let mut args = std::env::args();
    if args.len() == 3 {
        let folder = args.nth(1).unwrap();
        let ontology = Ontology::from_standard(&folder).unwrap();
        println!(
            "Ontology [{}] with {} terms",
            ontology.hpo_version(),
            ontology.len()
        );

        let filename = args.next().unwrap();
        let mut fh = File::create(filename).expect("Cannot create file");
        match fh.write_all(&ontology.as_bytes()) {
            Ok(_) => println!("Saved output"),
            Err(err) => println!("Error: {}", err),
        };
    } else {
        println!("Create a binary ontology from the JAX master data\n\n");
        println!("Usage\nobo_to_bin <SOURCE FOLDER> <OUTPUT FILENAME>");
        println!("e.g.:\nobo_to_bin example_data/ my_ontology.bin\n");
    }
}
