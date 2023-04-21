use std::path::Path;

use rayon::prelude::*;

use hpo::{HpoTermId, Ontology};

const N_INNER_ITEMS: usize = 10000;

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

/// Returns the pair of HpoTerms that contain the most common ancestors
fn multi_threaded_ancestors(ontology: &Ontology) -> (usize, HpoTermId, HpoTermId) {
    ontology
        .hpos()
        .par_bridge()
        .fold(
            || (0, HpoTermId::from(1u32), HpoTermId::from(1u32)),
            |outer_res, term1| {
                let inner_res = ontology
                    .hpos()
                    .take(N_INNER_ITEMS)
                    .filter(|term2| &term1 != term2)
                    .fold(
                        (0, HpoTermId::from(1u32), HpoTermId::from(1u32)),
                        |res, term2| {
                            let overlap = term1.common_ancestor_ids(&term2).len();
                            if overlap > res.0 {
                                (overlap, term1.id(), term2.id())
                            } else {
                                res
                            }
                        },
                    );
                if inner_res.0 > outer_res.0 {
                    inner_res
                } else {
                    outer_res
                }
            },
        )
        // rayon fold implementation requires another reduce call
        // to get to a single result
        .reduce(
            || (0, HpoTermId::from(1u32), HpoTermId::from(1u32)),
            |max, current| {
                if current.0 > max.0 {
                    current
                } else {
                    max
                }
            },
        )
}

/// Returns the pair of HpoTerms that contain the most common ancestors
fn single_treaded_ancestors(ontology: &Ontology) -> (usize, HpoTermId, HpoTermId) {
    ontology.hpos().fold(
        (0, HpoTermId::from(1u32), HpoTermId::from(1u32)),
        |outer_res, term1| {
            let inner_res = ontology
                .hpos()
                .take(N_INNER_ITEMS)
                .filter(|term2| &term1 != term2)
                .fold(
                    (0, HpoTermId::from(1u32), HpoTermId::from(1u32)),
                    |res, term2| {
                        let overlap = term1.common_ancestor_ids(&term2).len();
                        if overlap > res.0 {
                            (overlap, term1.id(), term2.id())
                        } else {
                            res
                        }
                    },
                );
            if inner_res.0 > outer_res.0 {
                inner_res
            } else {
                outer_res
            }
        },
    )
}

fn main() {
    let mut args = std::env::args();
    if args.len() < 2 {
        panic!("Usage: ./bench_read_ontology /path/to/ontology");
    }
    let ontology = read_ontology(&args.nth(1).unwrap());

    let res = if args.next().is_some() {
        multi_threaded_ancestors(&ontology)
    } else {
        single_treaded_ancestors(&ontology)
    };

    println!("Terms {} and {} have {} overlaps", res.1, res.2, res.0);
}
