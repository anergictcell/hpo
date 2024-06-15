//! Calculate the enrichment of diseases in a selected set of genes
//! and count how many genes are enriched with a p-value < 0.005 and > 2fold enrichment

use std::path::Path;

use rayon::prelude::*;

use hpo::stats::hypergeom::omim_disease_enrichment;
use hpo::Ontology;

const GENES: [&str; 6] = ["GBA1", "NPC1", "EZH2", "DMD", "MUC7", "ARID1B"];

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

fn multi_threaded_enrichment(ontology: &Ontology) -> usize {
    GENES
        .iter()
        .par_bridge()
        .fold(
            || 0,
            |total, gene| {
                let set = ontology
                    .gene_by_name(gene)
                    .unwrap_or_else(|| panic!("{gene} does not exist"))
                    .to_hpo_set(ontology);
                total
                    + omim_disease_enrichment(ontology, &set)
                        .iter()
                        .fold(0, |count, enrichment| {
                            if enrichment.pvalue() < 0.0000005 && enrichment.enrichment() > 100.0 {
                                count + 1
                            } else {
                                count
                            }
                        })
            },
        )
        .sum::<usize>()
}

fn single_threaded_enrichment(ontology: &Ontology) -> usize {
    GENES.iter().fold(0, |total, gene| {
        let set = ontology
            .gene_by_name(gene)
            .unwrap_or_else(|| panic!("{gene} does not exist"))
            .to_hpo_set(ontology);
        total
            + omim_disease_enrichment(ontology, &set)
                .iter()
                .fold(0, |count, enrichment| {
                    if enrichment.pvalue() < 0.0000005 && enrichment.enrichment() > 100.0 {
                        count + 1
                    } else {
                        count
                    }
                })
    })
}

fn main() {
    let mut args = std::env::args();
    if args.len() < 2 {
        panic!("Usage: ./bench_disease_enrichment /path/to/ontology");
    }
    let ontology = read_ontology(&args.nth(1).unwrap());

    let max = if args.next().is_some() {
        multi_threaded_enrichment(&ontology)
    } else {
        single_threaded_enrichment(&ontology)
    };

    println!("Highly enriched: {max}");
}
