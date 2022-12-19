use hpo::HpoTerm;
use hpo::HpoTermId;
use hpo::Ontology;

fn all_distances(ontology: &Ontology) {
    for term1 in ontology.hpos() {
        for term2 in ontology.hpos() {
            print_distance(&term1, &term2)
        }
    }
}

fn print_distance(term1: &HpoTerm, term2: &HpoTerm) {
    if let Some(dist) = term1.distance_to_ancestor(term2) {
        println!("{}\t{}\t{}", term1.id(), term2.id(), dist);
    }
}

fn main() {
    let ontology = Ontology::from_standard("./example_data/").unwrap();

    let mut args = std::env::args();
    if args.len() == 3 {
        let termid1 = HpoTermId::from(args.nth(1).unwrap());
        let termid2 = HpoTermId::from(args.next().unwrap());
        let term1 = ontology.hpo(&termid1).unwrap();
        let term2 = ontology.hpo(&termid2).unwrap();
        print_distance(&term1, &term2);
        if let Some(path) = term1.path_to_ancestor(&term2) {
            for term in path {
                println!("{}", term);
            }
        } else {
            println!("No common ancestor");
        }
    } else {
        all_distances(&ontology);
    }
}
