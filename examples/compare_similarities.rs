//! Compare 2 different similarity algorithms side by side and print
//! the result as a tab separated output.
//!
//! # Arguments:
//!
//! - Algorithm 1
//! - Algorithm 2
//! - Number of Hpo Terms to compare (leave blank to compare all)
//!
//! # Usage:
//! compare_similarities graphic mutation 100
//!

use std::{
    io::{self, Write},
    process,
};

use rayon::prelude::*;

use hpo::{similarity::Builtins, term::InformationContentKind, Ontology};
fn main() {
    let mut args = std::env::args();
    if args.len() < 3 {
        println!("Compare two different similarity algorithms\n\n");
        println!("Usage\ncompare_similarities <ALGORITHM> <ALGORITHM> [<NUMBER_OF_TERMS>]");
        println!("e.g.:\ncompare_similarities graphic resnik 100\n");
        process::exit(1)
    }
    let sim1_name = args.nth(1).expect("no algorithm 1 specified");
    let sim2_name = args.next().expect("no algorithm 2 specified");
    let n_comparisons = args
        .next()
        .map(|f| f.parse::<usize>().expect("number of terms must be numeric"))
        .unwrap_or(20_000);

    let ontology = Ontology::from_binary("./tests/ontology.hpo").unwrap();

    let ic_kind = InformationContentKind::Omim;
    let sim1 = Builtins::new(&sim1_name, ic_kind).expect("invalid algoritm 1");
    let sim2 = Builtins::new(&sim2_name, ic_kind).expect("invalid algoritm 2");

    let mut n_terms = 100_000;

    if let Some(n) = args.next() {
        if let Ok(items) = n.parse::<usize>() {
            n_terms = items;
        }
    }

    let scores: Vec<String> = ontology
        .into_iter()
        .take(n_terms)
        .par_bridge()
        .map(|term1| {
            let mut inner_score = Vec::new();
            for term2 in ontology.into_iter().take(n_comparisons) {
                let score1 = term1.similarity_score(&term2, &sim1);
                let score2 = term1.similarity_score(&term2, &sim2);
                inner_score.push(format!(
                    "{}\t{}\t{:.4}\t{:.4}\n",
                    term1.id(),
                    term2.id(),
                    score1,
                    score2
                ));
            }
            inner_score
        })
        .flatten()
        .collect();

    let mut stdout = io::stdout().lock();
    for row in scores {
        stdout.write_all(row.as_bytes()).unwrap();
    }
}
