use std::fs::File;
use std::io::BufRead;
use std::io::BufReader;
use std::time::SystemTime;

use hpo::parser;
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

    parser::phenotype_to_genes::parse("phenotype_to_genes.txt", collection);
    println!("finished linking genes");

    parser::phenotype_to_genes::parse("phenotype_to_genes.txt", collection);
    parser::phenotype_to_genes::parse("phenotype_to_genes.txt", collection);
    parser::phenotype_hpoa::parse("phenotype.hpoa", collection);
    collection.calculate_information_content();
    println!("finished IC calculation");
}

fn bench(collection: &Ontology, times: usize) {
    let mut count = 0;
    let mut terms: (HpoTermId, HpoTermId) =
        (HpoTermId::from("HP:0000001"), HpoTermId::from("HP:0000001"));
    let start = SystemTime::now();
    for term1 in collection.iter_terms() {
        for term2 in collection.iter_terms().take(times) {
            let overlap = term1.common_ancestor_ids(&term2).len();
            if overlap > count {
                count = overlap;
                terms = (*term1.id(), *term2.id());
            }
        }
    }
    let end = SystemTime::now();
    let duration = end.duration_since(start).unwrap();
    println!(
        "It took {} seconds for {} terms. {} and {} have {} overlaps",
        duration.as_secs(),
        std::cmp::min(times, collection.len()),
        terms.0,
        terms.1,
        count
    );
}

fn common_ancestors(termid1: HpoTermId, termid2: HpoTermId, ontology: &Ontology) {
    let term1 = ontology.get_term(&termid1).unwrap();
    let term2 = ontology.get_term(&termid2).unwrap();
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
    let mut collection = Ontology::default();
    from_file(&mut collection);

    println!("finished creating Ontology");

    let mut args = std::env::args();
    if args.len() == 3 {
       let term_id1 = args.nth(1).unwrap();
       let term_id2 = args.next().unwrap();
       common_ancestors(term_id1.into(), term_id2.into(), &collection);
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
