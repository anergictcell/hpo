use std::fs::File;
use std::io::BufRead;
use std::io::BufReader;
use std::time::SystemTime;

use hpo::GraphIc;
use hpo::parser;
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
    parser::phenotype_hpoa::parse("phenotype.hpoa", collection);
    collection.calculate_information_content();
    println!("Finished Information content calculation");
}

fn bench(collection: &Ontology, times: usize) {
    let start = SystemTime::now();
    let ic = GraphIc::new(hpo::InformationContentKind::Omim);
    for term1 in collection.iter_terms() {
        for term2 in collection.iter_terms().take(times) {
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
        collection.len(),
        std::cmp::min(times, collection.len()),
    );
}

fn main() {
    let mut collection = Ontology::default();
    from_file(&mut collection);

    println!("finished creating Ontology");

    bench(&collection, 5);
    bench(&collection, 50);
    bench(&collection, 500);
    bench(&collection, 1000);
    bench(&collection, 10000);
    bench(&collection, 20000);

    /*
    Expected times:
    It took 0 seconds for 17059 x 5 terms.
    It took 0 seconds for 17059 x 50 terms.
    It took 1 seconds for 17059 x 500 terms.
    It took 3 seconds for 17059 x 1000 terms.
    It took 42 seconds for 17059 x 10000 terms.
    It took 70 seconds for 17059 x 17059 terms.
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
