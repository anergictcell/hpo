//! Compare two different ontologies
//!
//! This is helpful for checking correctness of the parser modules
//! or to see changes after a new HPO release

use hpo::comparison::Comparison;
use hpo::{HpoTermId, Ontology};
use std::{path::Path, process};

fn ontology(path_arg: &str) -> Ontology {
    let path = Path::new(path_arg);

    match path.is_file() {
        true => Ontology::from_binary(path).unwrap(),
        false => Ontology::from_standard(&path.to_string_lossy()).unwrap(),
    }
}

/// Prints some basic stats about the differences
/// between two Ontologies
fn overview(diffs: &Comparison) {
    println!("#Numbers\n{}", diffs);
    println!("#Change\tID\tName");
    for term in diffs.added_hpo_terms() {
        println!("Added\t{}\t{}", term.id(), term.name());
    }
    for term in diffs.removed_hpo_terms() {
        println!("Removed\t{}\t{}", term.id(), term.name());
    }

    for gene in diffs.added_genes() {
        println!("Added\t{}\t{}", gene.id(), gene.name());
    }

    for gene in diffs.removed_genes() {
        println!("Removed\t{}\t{}", gene.id(), gene.name());
    }

    for disease in diffs.added_omim_diseases() {
        println!("Added\t{}\t{}", disease.id(), disease.name());
    }

    for disease in diffs.removed_omim_diseases() {
        println!("Removed\t{}\t{}", disease.id(), disease.name());
    }
}

/// Prints info about Term-specific changes
fn changed_terms(diffs: &Comparison) {
    println!(
        "#Term Delta\tID\tOld Name:New Name\tAdded Parents\tRemoved Parents\tObsolete\tReplacement"
    );
    for term in diffs.changed_hpo_terms() {
        print!("Delta\t{}", term.id());
        if let Some(names) = term.changed_name() {
            print!("\t{}:{}", names.0, names.1);
        } else {
            print!("\t.");
        }
        if let Some(added) = term.added_parents() {
            print!(
                "\t{}",
                added
                    .iter()
                    .map(|tid| tid.to_string())
                    .collect::<Vec<String>>()
                    .join(",")
            );
        } else {
            print!("\t.");
        }
        if let Some(removed) = term.removed_parents() {
            print!(
                "\t{}",
                removed
                    .iter()
                    .map(|tid| tid.to_string())
                    .collect::<Vec<String>>()
                    .join(",")
            );
        } else {
            print!("\t.");
        }
        print_obsolete_diff(term.changed_obsolete());
        print_replacement_diff(term.changed_replacement());
        println!();
    }
}

/// Prints info about Gene-specific changes
fn changed_genes(diffs: &Comparison) {
    println!("#Gene Delta\tID\tOld Name:New Name\tAdded Parents\tRemoved Parents");
    for term in diffs.changed_genes() {
        print!("Delta\t{}", term.id());
        if let Some(names) = term.changed_name() {
            print!("\t{}:{}", names.0, names.1);
        } else {
            print!("\t.");
        }
        if let Some(added) = term.added_terms() {
            print!(
                "\t{}",
                added
                    .iter()
                    .map(|tid| tid.to_string())
                    .collect::<Vec<String>>()
                    .join(",")
            );
        } else {
            print!("\t.");
        }
        if let Some(removed) = term.removed_terms() {
            print!(
                "\t{}",
                removed
                    .iter()
                    .map(|tid| tid.to_string())
                    .collect::<Vec<String>>()
                    .join(",")
            );
        } else {
            print!("\t.");
        }
        println!();
    }
}

/// Prints info about Gene-specific changes
fn changed_diseases(diffs: &Comparison) {
    println!("#Disease Delta\tID\tOld Name:New Name\tAdded Parents\tRemoved Parents");
    for term in diffs.changed_omim_diseases() {
        print!("Delta\t{}", term.id());
        if let Some(names) = term.changed_name() {
            print!("\t{}:{}", names.0, names.1);
        } else {
            print!("\t.");
        }
        if let Some(added) = term.added_terms() {
            print!(
                "\t{}",
                added
                    .iter()
                    .map(|tid| tid.to_string())
                    .collect::<Vec<String>>()
                    .join(",")
            );
        } else {
            print!("\t.");
        }
        if let Some(removed) = term.removed_terms() {
            print!(
                "\t{}",
                removed
                    .iter()
                    .map(|tid| tid.to_string())
                    .collect::<Vec<String>>()
                    .join(",")
            );
        } else {
            print!("\t.");
        }
        println!();
    }
}

fn print_obsolete_diff(obsoletes: Option<(bool, bool)>) {
    if let Some(obsolete) = obsoletes {
        print!(
            "\t{}:{}",
            if obsolete.0 { "obsolete" } else { "in-use" },
            if obsolete.1 { "obsolete" } else { "in-use" }
        );
    } else {
        print!("\t.");
    }
}

fn print_replacement_diff(replacements: Option<(Option<HpoTermId>, Option<HpoTermId>)>) {
    if let Some(replacements) = replacements {
        print!(
            "\t{}:{}",
            replacements
                .0
                .map(|id| id.to_string())
                .unwrap_or("-".to_string()),
            replacements
                .1
                .map(|id| id.to_string())
                .unwrap_or("-".to_string())
        );
    } else {
        print!("\t.");
    }
}

fn main() {
    let mut args = std::env::args();

    if args.len() != 3 {
        println!("Compare two Ontologies to each other and print the differences\n\n");
        println!("Usage:\ncompare_ontologies </PATH/TO/ONTOLOGY> </PATH/TO/OTHER-ONTOLOGY>");
        println!("e.g.:\ncompare_ontologies tests/ontology.hpo tests/ontology_v2.hpo:\n");
        process::exit(1)
    }
    let arg_old = args.nth(1).unwrap();
    let arg_new = args.next().unwrap();

    let lhs = ontology(&arg_old);
    let rhs = ontology(&arg_new);

    let diffs = lhs.compare(&rhs);

    overview(&diffs);
    changed_terms(&diffs);
    changed_genes(&diffs);
    changed_diseases(&diffs);
}
