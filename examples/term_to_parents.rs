use std::fs::File;
use std::io::BufRead;
use std::io::BufReader;


use hpo::Ontology;

fn from_file(collection: &mut Ontology) {
    let file = File::open("terms.txt").unwrap();
    let reader = BufReader::new(file);
    for term in reader.lines() {
        collection.add_term_by_name(&term.unwrap());
    }

    let file = File::open("connections.txt").unwrap();
    let reader = BufReader::new(file);
    for line in reader.lines() {
        let line = line.unwrap();
        let cols: Vec<&str> = line.splitn(2, '\t').collect();
        collection.add_parent(cols[1].into(), cols[0].into());
    }

    collection.create_cache();

    let file = File::open("phenotype_to_genes.txt").unwrap();
    let reader = BufReader::new(file);
    for line in reader.lines() {
        let line = line.unwrap();
        if line.starts_with('#') { continue; }
        let cols: Vec<&str> = line.trim().split('\t').collect();
        let gene_id = collection.add_gene(cols[3], cols[2]).unwrap();
        collection.link_gene_term(&cols[0].into(), gene_id);
    }
}

fn main() {
    let mut collection = Ontology::default();
    from_file(&mut collection);

    let mut parents: Vec<String> = Vec::new();
    for term in collection.iter_terms() {
        parents.clear();
        for parent in term.all_parents() {
            parents.push(parent.id().to_string());
        }
        parents.sort();
        println!("{}\t{}", term.id(), parents.join(","));
    }
}

/*
Python comparison code
from pyhpo import Ontology

_ = Ontology()

with open("term2parents.py.txt", "w") as fh:
    for term in Ontology:
        parents = []
        for parent in term.all_parents:
            parents.append(parent.id)
        parents = ",".join(sorted(parents))
        _ = fh.write(f"{term.id}\t{parents}\n")

*/