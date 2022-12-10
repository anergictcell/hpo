use hpo::HpoTerm;
use hpo::HpoTermId;
use hpo::Ontology;

fn bench(ontology: &Ontology) {
    for term1 in ontology.hpos() {
        for term2 in ontology.hpos() {
            print_distance(&term1, &term2)
        }
    }
}

fn print_distance(term1: &HpoTerm, term2: &HpoTerm) {
    if let Some(dist) = term1.distance_to_term(term2) {
        println!("{}\t{}\t{}", term1.id(), term2.id(), dist);
    }
}

fn main() {
    let ontology = Ontology::from_standard("./example_data/");

    println!("finished creating Ontology");

    let mut args = std::env::args();
    if args.len() == 3 {
        let termid1 = HpoTermId::from(args.nth(1).unwrap());
        let termid2 = HpoTermId::from(args.next().unwrap());
        let term1 = ontology.hpo(&termid1).unwrap();
        let term2 = ontology.hpo(&termid2).unwrap();
        print_distance(&term1, &term2);
        for t in term1.common_ancestors(&term2) {
            println!("{}", t.id());
        }
        println!("----");
        println!(">> {} <<", term1.id());
        for t in term1.all_parents() {
            println!("{}", t.id());
        }
        println!("----");
        println!(">> {} <<", term2.id());
        for t in term2.all_parents() {
            println!("{}", t.id());
        }
    } else {
        bench(&ontology);
    }
}
