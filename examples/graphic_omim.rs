use std::fs::File;
use std::io::BufRead;
use std::io::BufReader;

use hpo::annotations::OmimDiseaseId;
use hpo::parser;
use hpo::term::HpoTermIterator;
use hpo::GraphIc;
use hpo::HpoTerm;
use hpo::Ontology;

fn from_file(collection: &mut Ontology) {
    let file = File::open("terms.txt").unwrap();
    let reader = BufReader::new(file);
    for term in reader.lines() {
        collection.add_term_by_name(&term.unwrap());
    }
    collection.shrink_to_fit();

    let file = File::open("connections.txt").unwrap();
    let reader = BufReader::new(file);
    for line in reader.lines() {
        let line = line.unwrap();
        let cols: Vec<&str> = line.splitn(2, '\t').collect();
        collection.add_parent(cols[1].into(), cols[0].into());
    }
    collection.create_cache();

    parser::phenotype_to_genes::parse("phenotype_to_genes.txt", collection);
    parser::phenotype_to_genes::parse("phenotype_to_genes.txt", collection);
    parser::phenotype_hpoa::parse("phenotype.hpoa", collection);
    collection.calculate_information_content();
}

fn scores_for_disease(omimid: OmimDiseaseId, ontology: &Ontology) {
    let ic = GraphIc::new(hpo::InformationContentKind::Omim);

    let disease = ontology.omim_disease(&omimid).unwrap();
    println!("Using Disease [{}] {}", disease.id(), disease.name());
    let terms: Vec<HpoTerm> = HpoTermIterator::new(disease.hpo_terms(), ontology).collect();
    for term1 in &terms {
        for term2 in &terms {
            let overlap = term1.similarity_score(term2, &ic);
            println!(
                "{}\t{}\t{}",
                term1.id(),
                term2.id(),
                (overlap * 10_000_000.0).floor()
            );
        }
    }
}

fn main() {
    let mut collection = Ontology::default();
    from_file(&mut collection);

    scores_for_disease("300486".try_into().unwrap(), &collection);

    /*
    Expected times:
    It took 0 seconds for 5 terms. HP:0000001 and HP:0000001 have 1 score
    It took 0 seconds for 50 terms. HP:0000001 and HP:0000001 have 1 score
    It took 1 seconds for 500 terms. HP:0000001 and HP:0000001 have 1 score
    It took 3 seconds for 1000 terms. HP:0000001 and HP:0000001 have 1 score
    It took 40 seconds for 10000 terms. HP:0000001 and HP:0000001 have 1 score
    */
}

/*
Python comparison code
from pyhpo import Ontology
from pyhpo.annotations import Omim

_ = Ontology()

disease = Omim.get(300486)
# OmimDisease(id=300486, name='Mental retardation, X-linked, with cerebellar hypoplasia and distinctive facial appearance', hpo={256, 3593, 1290, 1419, 400, 2066, 276, 7065, 28, 6951, 1320, 1321, 46, 303, 30260, 54, 1344, 31936, 448, 322, 2119, 336, 11220, 2007, 601, 219, 733, 1249, 1250, 486, 742, 2280, 744, 490, 1257, 1263, 752, 3189, 25336, 639}, negative_hpo=set(), diseasetype='Omim')

terms = [Ontology[idx] for idx in disease.hpo]

with open("graphic_omim.py.txt", "w") as fh:
    _ = fh.write(f"Using Disease [{disease.id}] {disease.name}\n")
    for term1 in terms:
        for term2 in terms:
            ic = term1.similarity_score(term2, kind="omim", method="graphic")
            _ = fh.write(f"{term1.id}\t{term2.id}\t{int(ic * 10_000_000)}\n")

*/
