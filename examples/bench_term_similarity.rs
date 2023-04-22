use std::path::Path;

use rayon::prelude::*;

use hpo::{
    similarity::{Builtins, Similarity},
    Ontology,
};

const N_COMPARISONS: usize = 1000;

/// Constructs an Ontology from either JAX-supplied obo data
/// or from an HPO-specific binary file
fn read_ontology(path_arg: &str) -> Ontology {
    let path = Path::new(path_arg);

    if path.is_file() {
        Ontology::from_binary(path).unwrap()
    } else {
        Ontology::from_standard(&path.to_string_lossy()).unwrap()
    }
}

fn multi_threaded_similarity(ontology: &Ontology) -> i32 {
    let sim = Builtins::GraphIc(hpo::term::InformationContentKind::Omim);
    ontology
        .hpos()
        .par_bridge()
        .map(|term_1| {
            ontology
                .hpos()
                .take(N_COMPARISONS)
                .fold(0, |count, term_2| {
                    if sim.calculate(&term_1, &term_2) > 0.9 {
                        count + 1
                    } else {
                        count
                    }
                })
        })
        .sum()
}

fn single_treaded_similarity(ontology: &Ontology) -> i32 {
    let sim = Builtins::GraphIc(hpo::term::InformationContentKind::Omim);
    let mut count = 0;
    for term_1 in ontology.hpos() {
        for term_2 in ontology.hpos().take(N_COMPARISONS) {
            if sim.calculate(&term_1, &term_2) > 0.9 {
                count += 1;
            }
        }
    }
    count
}

fn main() {
    let mut args = std::env::args();
    if args.len() < 2 {
        panic!("Usage: ./bench_read_ontology /path/to/ontology");
    }
    let ontology = read_ontology(&args.nth(1).unwrap());

    let count = if args.next().is_some() {
        multi_threaded_similarity(&ontology)
    } else {
        single_treaded_similarity(&ontology)
    };

    println!("Comparisons above 0.9 {count}");
}
