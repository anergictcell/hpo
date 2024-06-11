use hpo::{HpoTermId, Ontology};
use std::io::Write;
use std::process;
use std::{fs::File, path::Path};

fn ontology(path_arg: &str) -> Ontology {
    let path = Path::new(path_arg);

    match path.is_file() {
        true => Ontology::from_binary(path).unwrap(),
        false => Ontology::from_standard(&path.to_string_lossy()).unwrap(),
    }
}

fn main() {
    let mut args = std::env::args();

    if args.len() != 3 {
        println!("Create a subontology\n\n");
        println!("Usage:\nsubontology </PATH/TO/ONTOLOGY> <COMMA SEPARATED LIST OF TERM IDs>");
        println!("e.g.:\nsubontology tests/ontology.hpo HP:0000001,HP:0000005,HP:0000007,HP:0000118,HP:0000707,HP:0000818,HP:0000864,HP:0001939,HP:0002011,HP:0002715,HP:0003581,HP:0003674,HP:0010662,HP:0010978,HP:0011017,HP:0012285,HP:0012443,HP:0012638,HP:0012639,HP:0012647,HP:0012648,HP:0012823,HP:0025454,HP:0031797,HP:0034345,HP:0100547\n");
        process::exit(1)
    }
    let path = args.nth(1).unwrap();
    let terms = args.next().unwrap();
    let terms = terms
        .split(',')
        .map(|id| HpoTermId::try_from(id).expect("Term-ID must be properly formatted"));

    let ontology = ontology(&path);

    let sub = ontology
        .sub_ontology(
            ontology
                .hpo(1u32)
                .expect("Root term must exist in ontology"),
            terms.map(|id| ontology.hpo(id).expect("Term must exist in ontology")),
        )
        .expect("sub ontology must work");

    let mut fh = File::create("out.hpo").expect("Cannot create file");
    match fh.write_all(&sub.as_bytes()) {
        Ok(_) => println!("Saved subontology to out.hpo"),
        Err(err) => println!("Error: {}", err),
    };
}
