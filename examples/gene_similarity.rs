use rayon::prelude::*;
use std::io;
use std::io::Write;
use std::{env::Args, time::SystemTime};

use hpo::{
    annotations::{Gene, GeneId},
    similarity::{GraphIc, GroupSimilarity, StandardCombiner},
    term::HpoGroup,
    HpoSet, HpoTermId, Ontology,
};

/// Calculates the similarity score of two diseases
/// The two diseases are specified as OMIM-ID via CLI arguments
fn compare_two_genes(
    ontology: &Ontology,
    sim: &GroupSimilarity<GraphIc, StandardCombiner>,
    mut args: Args,
) {
    let symbol_a = args.nth(1).unwrap();
    let symbol_b = args.next().unwrap();

    let gene_a = ontology
        .gene_by_name(&symbol_a)
        .expect("The first gene is not part of the Ontology");
    let gene_b = ontology
        .gene_by_name(&symbol_b)
        .expect("The second gene is not part of the Ontology");

    let set_a = gene_a.to_hpo_set(ontology);
    let set_b = gene_b.to_hpo_set(ontology);

    let res = sim.calculate(&set_a, &set_b);
    println!("Similarity is {res}");
}

/// Calculates the similarity score of a custom HPO-Set
/// to all OMIM diseases
/// The HPO-Set is specified as a comma separated list of HPO-Term-IDs
fn compare_custom_set_to_genes(
    ontology: &Ontology,
    sim: &GroupSimilarity<GraphIc, StandardCombiner>,
    terms: String,
) {
    let hpo_terms = terms.split(',');
    let mut group = HpoGroup::default();
    for t in hpo_terms {
        group.insert(HpoTermId::try_from(t).expect("Invalid HpoTermId"));
    }
    let set_a = HpoSet::new(ontology, group);

    let start = SystemTime::now();
    let mut results: Vec<(&Gene, f32)> = ontology
        .genes()
        .par_bridge()
        .map(|gene| {
            let res = sim.calculate(&set_a, &gene.to_hpo_set(ontology));
            (gene, res)
        })
        .collect();
    let end = SystemTime::now();
    let duration = end.duration_since(start).unwrap();

    println!(
        "Number of comparisons: {} in {} sec",
        results.len(),
        duration.as_secs()
    );
    results.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());

    let mut stdout = io::stdout().lock();
    for x in results {
        stdout
            .write_all(format!("{}\t{}\t{}\n", x.0.id(), x.0.name(), x.1).as_bytes())
            .unwrap();
    }
}

/// Calculate the pairwise similarity of all diseases to <num> other
/// diseases.
fn cross_compare_genes(
    ontology: &Ontology,
    sim: &GroupSimilarity<GraphIc, StandardCombiner>,
    num: usize,
) {
    let start = SystemTime::now();
    let results: Vec<(GeneId, GeneId, f32)> = ontology
        .genes()
        .par_bridge()
        .flat_map(|gene_a| {
            ontology
                .genes()
                .take(num)
                .map(|gene_b| {
                    let res =
                        sim.calculate(&gene_a.to_hpo_set(ontology), &gene_b.to_hpo_set(ontology));
                    (*gene_a.id(), *gene_b.id(), res)
                })
                .collect::<Vec<(GeneId, GeneId, f32)>>()
        })
        .collect();
    let end = SystemTime::now();
    let duration = end.duration_since(start).unwrap();

    println!(
        "Number of comparisons: {} in {} sec",
        results.len(),
        duration.as_secs()
    );

    let mut stdout = io::stdout().lock();
    for x in results {
        stdout
            .write_all(format!("{}\t{}\t{}\n", x.0, x.1, x.2).as_bytes())
            .unwrap();
    }
}

fn main() {
    let ontology = Ontology::from_binary("tests/example.hpo").unwrap();
    let sim = GroupSimilarity::default();

    let mut args = std::env::args();

    match args.len() {
        3 => compare_two_genes(&ontology, &sim, args),
        2 => {
            let arg = args.nth(1).unwrap();
            if let Ok(num) = arg.parse::<usize>() {
                // integer provided, using disease x disease comparisons
                cross_compare_genes(&ontology, &sim, num);
            } else {
                // List of HPO terms provided
                compare_custom_set_to_genes(&ontology, &sim, arg);
            }
        }
        _ => {
            println!("Calculate similarities of genes\n\n");
            println!("\
                There are 3 different options:\n\
                - Compare 2 genes\n\
                    gene_similarity <Gene symbol> <Gene symbol>\n\
                    gene_similarity 618395 615368\n\n\
                - Compare all genes to a custom HPO-Set\n\
                    gene_similarity <HPO-TERM-IDs>\n\
                    gene_similarity HP:0000750,HP:0000752,HP:0001249,HP:0007018,HP:0010818,HP:0011463\n\n\
                - Cross compare N genes\n\
                    gene_similarity <NUMBER OF COMPARISONS>\n\
                    gene_similarity 20
            ");
        }
    }
}
