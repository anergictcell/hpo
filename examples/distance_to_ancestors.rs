use std::fs::File;
use std::io::BufRead;
use std::io::BufReader;
use hpo::HpoTerm;
use hpo::HpoTermId;
use hpo::Ontology;

fn from_file(collection: &mut Ontology) {
    let file = File::open("terms.txt").unwrap();
    let reader = BufReader::new(file);
    for term in reader.lines() {
        collection.add_term_by_name(&term.unwrap());
    }
    collection.shrink_to_fit();
    let file = File::open("connections.txt").unwrap();
    let reader = BufReader::new(file);
    for line in reader.lines() {
        let line = line.unwrap();
        let cols: Vec<&str> = line.splitn(2, '\t').collect();
        collection.add_parent(cols[1].into(), cols[0].into());
    }
    collection.create_cache();
}

fn bench(collection: &Ontology, times: usize) {
    for term1 in collection.iter_terms() {
        for term2 in collection.iter_terms().take(times) {
            print_distance(&term1, &term2)
        }
    }
}

fn print_distance(term1: &HpoTerm, term2: &HpoTerm) {
    let dist = term1.distance_to_ancestor(term2).unwrap_or(1_000_000);
    println!("{}\t{}\t{}", term1.id(), term2.id(), dist);
}

fn main() {
    let mut collection = Ontology::default();
    from_file(&mut collection);

    println!("finished creating Ontology");

    let mut args = std::env::args();
    if args.len() == 3 {
        let termid1 = HpoTermId::from(args.nth(1).unwrap());
        let termid2 = HpoTermId::from(args.next().unwrap());
        let term1 = collection.get_term(&termid1).unwrap();
        let term2 = collection.get_term(&termid2).unwrap();
        print_distance(&term1, &term2);
        if let Some(path) = term1.path_to_ancestor(&term2) {
            for term in path {
                println!("{}", term);
            }
        } else {
            println!("No common ancestor");
        }
    } else {
        bench(&collection, 500);
        bench(&collection, 1000);
        bench(&collection, 10000);
        bench(&collection, 20000);
    }


    /*
    Expected times:
    It took 0 seconds for 500 terms. HP:0007768 and HP:0000631 have 24 overlaps
    It took 1 seconds for 1000 terms. HP:0005617 and HP:0001215 have 35 overlaps
    It took 13 seconds for 10000 terms. HP:0009640 and HP:0009640 have 42 overlaps
    It took 22 seconds for 17059 terms. HP:0009640 and HP:0009640 have 42 overlaps
    */
}

/*
Python comparison code
from pyhpo import Ontology
import time

_ = Ontology()

time1 = time.time()
common = 0
terms = (None, None)
for term1 in Ontology:
    # for term2 in list(Ontology)[0:5000]:
    for term2 in Ontology:
            overlap = term1.common_ancestors(term2)
            if len(overlap) > common:
                    common = len(overlap)
                    terms = (term1, term2)

time2 = time.time()
print('It took {:.3f} s. {} and {} have {} overlaps '.format(
    time2-time1,
    terms[0].id,
    terms[1].id,
    common
))

# 5000:
# It took 82.564 s. HP:0004223 and HP:0004223 have 41 overlaps

# All terms:
# It took 298.564 s. HP:0009640 and HP:0009640 have 42 overlaps

*/
