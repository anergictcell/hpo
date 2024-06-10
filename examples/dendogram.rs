use std::env::Args;
use std::path::Path;

use hpo::annotations::{Disease, OmimDisease};
use rayon::prelude::*;

use hpo::similarity::GroupSimilarity;
use hpo::stats::Linkage;
use hpo::utils::Combinations;
use hpo::HpoSet;
use hpo::Ontology;

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

fn single_threaded_distance(combs: Combinations<HpoSet<'_>>) -> Vec<f32> {
    let sim = GroupSimilarity::default();
    combs
        .map(|comp| 1.0 - sim.calculate(comp.0, comp.1))
        .collect()
}

fn multi_threaded_distance(combs: Combinations<HpoSet<'_>>) -> Vec<f32> {
    let sim = GroupSimilarity::default();
    let x: Vec<(&HpoSet, &HpoSet)> = combs.collect();
    x.par_iter()
        .map(|comp| 1.0 - sim.calculate(comp.0, comp.1))
        .collect()
}

fn dendogram<'a>(ontology: &'a Ontology, args: &mut Args) -> (Linkage<'a>, Vec<&'a str>) {
    let diseases: Vec<&OmimDisease> = ontology
        .omim_diseases()
        .take(args.next().unwrap().parse::<usize>().unwrap())
        .collect();

    let sets = diseases
        .iter()
        .map(|g| g.to_hpo_set(ontology).without_modifier().child_nodes());
    let linkage = if args.next().is_some() {
        Linkage::single(sets, multi_threaded_distance)
    } else {
        Linkage::single(sets, single_threaded_distance)
    };

    let sorted_diseases = linkage
        .indicies()
        .iter()
        .map(|idx| diseases.get(*idx).unwrap().name())
        .collect::<Vec<&str>>();
    (linkage, sorted_diseases)
}

fn main() {
    let mut args = std::env::args();
    if args.len() < 3 {
        panic!("Usage: ./dendogram /path/to/ontology NUM_DISEASES");
    }
    let ontology = read_ontology(&args.nth(1).unwrap());

    let res = dendogram(&ontology, &mut args);
    for cluster in res.0.cluster() {
        println!("{:?}", cluster);
    }
    for d in res.1 {
        println!("{d}");
    }
}
