use std::path::Path;

use hpo::Ontology;

fn main() {
    let mut args = std::env::args();
    if args.len() < 2 {
        panic!("Usage: ./bench_read_ontology /path/to/ontology");
    }
    let path_arg = args.nth(1).unwrap();
    let path = Path::new(&path_arg);

    let ontology = match path.is_file() {
        true => Ontology::from_binary(path).unwrap(),
        false => Ontology::from_standard(&path.to_string_lossy()).unwrap(),
    };
    println!("Loaded the Ontology with {} terms", ontology.len());
}
