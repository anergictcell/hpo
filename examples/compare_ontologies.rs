use hpo::{HpoTermId, Ontology};
use std::path::Path;

fn main() {
    let mut args = std::env::args();

    if args.len() != 3 {
        panic!("Usage: ./compare_ontologies /path/to/old/ontology /path/to/new/ontology");
    }
    let arg_old = args.nth(1).unwrap();
    let arg_new = args.next().unwrap();

    let path_old = Path::new(&arg_old);
    let path_new = Path::new(&arg_new);

    let lhs = match path_old.is_file() {
        true => Ontology::from_binary(path_old).unwrap(),
        false => Ontology::from_standard(&path_old.to_string_lossy()).unwrap(),
    };

    let rhs = match path_new.is_file() {
        true => Ontology::from_binary(path_new).unwrap(),
        false => Ontology::from_standard(&path_new.to_string_lossy()).unwrap(),
    };

    let diffs = lhs.compare(&rhs);
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
