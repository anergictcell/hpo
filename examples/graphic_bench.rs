use std::env::args;
use std::time::SystemTime;

use hpo::HpoTerm;
use hpo::HpoTermId;
use rayon::prelude::*;

use hpo::GraphIc;
use hpo::Ontology;

fn bench(ontology: &Ontology, times: usize) {
    let start = SystemTime::now();
    let ic = GraphIc::new(hpo::InformationContentKind::Omim);
    for term1 in ontology.hpos() {
        for term2 in ontology.hpos().take(times) {
            let overlap = term1.similarity_score(&term2, &ic);
            if overlap > 1.1 {
                println!("This part is never reached but is left so that the compiler doesn't optimize the loop away :)")
            }
        }
    }
    let end = SystemTime::now();
    let duration = end.duration_since(start).unwrap();
    println!(
        "It took {} seconds for {} x {} terms.",
        duration.as_secs(),
        ontology.len(),
        std::cmp::min(times, ontology.len()),
    );
}

fn parallel(ontology: &Ontology, times: usize) {
    let start = SystemTime::now();
    let ic = GraphIc::new(hpo::InformationContentKind::Omim);

    let terms: Vec<HpoTerm> = ontology.into_iter().collect();
    let scores: Vec<(HpoTermId, HpoTermId, f32)> = terms.par_iter()
    .map(|term1| {
        let mut inner_score = Vec::new();
        for term2 in ontology.into_iter().take(times) {
            let overlap = term1.similarity_score(&term2, &ic);
            inner_score.push((*term1.id(), *term2.id(), overlap));
            if overlap > 1.1 {
                println!("This part is never reached but is left so that the compiler doesn't optimize the loop away :)")
            }
        }
        inner_score
    }).flatten().collect();

    let end = SystemTime::now();
    let duration = end.duration_since(start).unwrap();
    println!(
        "It took {} seconds for {} x {} terms: {}.",
        duration.as_secs(),
        ontology.len(),
        std::cmp::min(times, ontology.len()),
        scores.len()
    );
}

fn main() {
    let ontology = Ontology::from_standard("./example_data/").unwrap();

    println!("finished creating Ontology");

    if args().len() == 2 {
        parallel(&ontology, 5);
        parallel(&ontology, 50);
        parallel(&ontology, 500);
        parallel(&ontology, 1000);
        parallel(&ontology, 10000);
        parallel(&ontology, 20000);
    } else {
        bench(&ontology, 5);
        bench(&ontology, 50);
        bench(&ontology, 500);
        bench(&ontology, 1000);
        bench(&ontology, 10000);
        bench(&ontology, 20000);
    }

    /*
    Expected times (single threaded):
    It took 0 seconds for 17059 x 5 terms.
    It took 0 seconds for 17059 x 50 terms.
    It took 1 seconds for 17059 x 500 terms.
    It took 3 seconds for 17059 x 1000 terms.
    It took 42 seconds for 17059 x 10000 terms.
    It took 70 seconds for 17059 x 17059 terms.


    Expected times (using rayon):
    It took 0 seconds for 17059 x 5 terms: 85295.
    It took 0 seconds for 17059 x 50 terms: 852950.
    It took 0 seconds for 17059 x 500 terms: 8529500.
    It took 0 seconds for 17059 x 1000 terms: 17059000.
    It took 8 seconds for 17059 x 10000 terms: 170590000.
    It took 13 seconds for 17059 x 17059 terms: 291009481.
    */
}

/*
Python comparison code
from pyhpo import Ontology
import time

_ = Ontology()

terms = 100

time1 = time.time()
for term1 in list(Ontology)[0:terms]:
    # for term2 in list(Ontology)[0:terms]:
    for term2 in Ontology:
        ic = term1.similarity_score(term2, kind="omim", method="graphic")
        if ic > 1.1:
            print("This should never happen")
            print("Term1: {term1}, 2: {term2}, Score: {}", ic)

time2 = time.time()
print('It took {:.3f} s for {} x {} terms.'.format(
    time2-time1,
    terms,
    terms
))

# 100:
# It took 9.542 s for 100 x 100 terms.

# 200:
# It took 19.583 s for 200 x 200 terms.

# 400
# It took 38.634 s for 400 x 400 terms.

*/
