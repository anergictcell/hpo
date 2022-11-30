use std::fs::File;
use std::io::BufRead;
use std::io::BufReader;
use std::time::SystemTime;

use hpo::HpoTermId;
use hpo::Ontology;

fn from_file(collection: &mut Ontology) {
    println!("Reading terms");
    let file = File::open("terms.txt").unwrap();
    let reader = BufReader::new(file);
    for term in reader.lines() {
        collection.add_term_by_name(&term.unwrap());
    }
    collection.shrink_to_fit();
    println!("finished adding terms");

    let file = File::open("connections.txt").unwrap();
    let reader = BufReader::new(file);
    for line in reader.lines() {
        let line = line.unwrap();
        let cols: Vec<&str> = line.splitn(2, '\t').collect();
        collection.add_parent(cols[1].into(), cols[0].into());
    }
    println!("finished adding connections");
    collection.create_cache();
    println!("finished caching");

    let file = File::open("phenotype_to_genes.txt").unwrap();
    let reader = BufReader::new(file);
    for line in reader.lines() {
        let line = line.unwrap();
        if line.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = line.trim().split('\t').collect();
        let gene_id = collection.add_gene(cols[3], cols[2]).unwrap();
        collection.link_gene_term(&cols[0].into(), gene_id);
    }
    println!("finished linking genes");
}

fn bench(collection: &Ontology, times: usize) {
    let mut count = 0;
    let mut terms: (HpoTermId, HpoTermId) =
        (HpoTermId::from("HP:0000001"), HpoTermId::from("HP:0000001"));
    let start = SystemTime::now();
    for term1 in collection.iter_terms() {
        for term2 in collection.iter_terms().take(times) {
            let overlap = term1.overlap(&term2).len();
            if overlap > count {
                count = overlap;
                terms = (*term1.id(), *term2.id());
            }
        }
    }
    let end = SystemTime::now();
    let duration = end.duration_since(start).unwrap();
    println!(
        "Term: it took {} seconds for {} terms. {} and {} have {} overlaps",
        duration.as_secs(),
        std::cmp::min(times, collection.len()),
        terms.0,
        terms.1,
        count
    );
}

fn main() {
    let mut collection = Ontology::default();
    from_file(&mut collection);

    println!("finished creating Ontology");

    bench(&collection, 500);
    bench(&collection, 1000);
    bench(&collection, 10000);
    bench(&collection, 20000);
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
