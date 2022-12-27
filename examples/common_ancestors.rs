use std::time::SystemTime;

use hpo::HpoTermId;
use hpo::Ontology;

fn bench(ontology: &Ontology, times: usize) {
    let mut count = 0;
    let mut terms: (HpoTermId, HpoTermId) = (
        HpoTermId::try_from("HP:0000001").unwrap(),
        HpoTermId::try_from("HP:0000001").unwrap(),
    );
    let start = SystemTime::now();
    ontology.hpos().for_each(|term1| {
        for term2 in ontology.hpos().take(times) {
            let overlap = term1.common_ancestor_ids(&term2).len();
            if overlap > count {
                count = overlap;
                terms = (*term1.id(), *term2.id());
            }
        }
    });
    let end = SystemTime::now();
    let duration = end.duration_since(start).unwrap();
    println!(
        "It took {} seconds for {} terms. {} and {} have {} overlaps",
        duration.as_secs(),
        std::cmp::min(times, ontology.len()),
        terms.0,
        terms.1,
        count
    );
}

fn common_ancestors(termid1: HpoTermId, termid2: HpoTermId, ontology: &Ontology) {
    let term1 = ontology.hpo(&termid1).unwrap();
    let term2 = ontology.hpo(&termid2).unwrap();
    for term in term1.common_ancestors(&term2) {
        println!(
            "Term {} | IC {} | nOmim {} | nGene {}",
            term.id(),
            term.information_content().omim_disease(),
            term.omim_diseases().count(),
            term.genes().count()
        );
    }
}

fn main() {
    let ontology = Ontology::from_standard("./example_data/").unwrap();

    let mut args = std::env::args();
    if args.len() == 3 {
        let term_id1 = args.nth(1).unwrap();
        let term_id2 = args.next().unwrap();
        common_ancestors(term_id1.into(), term_id2.into(), &ontology);
    } else {
        bench(&ontology, 500);
        bench(&ontology, 1000);
        bench(&ontology, 10000);
        bench(&ontology, 20000);
    }

    /*
    Expected times:
    It took 0 seconds for 500 terms. HP:0007768 and HP:0000631 have 23 overlaps
    It took 1 seconds for 1000 terms. HP:0030675 and HP:0001215 have 34 overlaps
    It took 12 seconds for 10000 terms. HP:0009640 and HP:0009640 have 41 overlaps
    It took 21 seconds for 17059 terms. HP:0009640 and HP:0009640 have 41 overlaps
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
    for term2 in list(Ontology)[0:1000]:
    #for term2 in Ontology:
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
