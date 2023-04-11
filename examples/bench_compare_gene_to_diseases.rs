use std::path::Path;

use rayon::prelude::*;

use hpo::similarity::GroupSimilarity;
use hpo::{HpoSet, Ontology};

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

fn multi_threaded_comparison(ontology: &Ontology, gene: &HpoSet) -> usize {
    let sim = GroupSimilarity::default();
    ontology
        .omim_diseases()
        .par_bridge()
        .map(|disease| {
            if sim.calculate(&disease.to_hpo_set(ontology), gene) > 0.6 {
                1
            } else {
                0
            }
        })
        .sum::<usize>()
}

fn single_threaded_comparison(ontology: &Ontology, gene: &HpoSet) -> usize {
    let sim = GroupSimilarity::default();
    let mut max = 0;
    for disease in ontology.omim_diseases() {
        if sim.calculate(&disease.to_hpo_set(ontology), gene) > 0.6 {
            max += 1;
        }
    }
    max
}

fn main() {
    let mut args = std::env::args();
    if args.len() < 2 {
        panic!("Usage: ./bench_sets_from_annotations /path/to/ontology");
    }
    let ontology = read_ontology(&args.nth(1).unwrap());

    let gba = ontology.gene_by_name("GBA1").unwrap().to_hpo_set(&ontology);

    let max = if args.next().is_some() {
        multi_threaded_comparison(&ontology, &gba)
    } else {
        single_threaded_comparison(&ontology, &gba)
    };

    println!("Highly similar: {max}");
}
