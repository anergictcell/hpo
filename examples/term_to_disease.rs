//! Prints every term and its associated Omim and Orpha Disease

use hpo::annotations::Disease;
use hpo::Ontology;

fn main() {
    let ontology = Ontology::from_binary("tests/ontology.hpo").unwrap();

    let mut diseases: Vec<String> = Vec::new();
    for term in ontology.hpos() {
        diseases.clear();
        for disease in term.omim_diseases() {
            diseases.push(disease.id().to_string());
        }
        diseases.sort();
        println!("{}\t{}", term.id(), diseases.join(","));
    }

    for term in ontology.hpos() {
        diseases.clear();
        for disease in term.orpha_diseases() {
            diseases.push(disease.id().to_string());
        }
        diseases.sort();
        println!("{}\t{}", term.id(), diseases.join(","));
    }
}

/*
Python comparison code
from pyhpo import Ontology

_ = Ontology()

with open("term2omim.py.txt", "w") as fh:
    for term in Ontology:
        omim_diseases = []
        for disease in term.omim_diseases:
            omim_diseases.append(f"OMIM:{disease.id}")
        omim_diseases = ",".join(sorted(omim_diseases))
        _ = fh.write(f"{term.id}\t{omim_diseases}\n")

    for term in Ontology:
        orpha_diseases = []
        for disease in term.orpha_diseases:
            orpha_diseases.append(f"ORPHA:{disease.id}")
        orpha_diseases = ",".join(sorted(orpha_diseases))
        _ = fh.write(f"{term.id}\t{orpha_diseases}\n")

*/

/*
diff-ing

diff \
<(sort example_data/term2omim.rs.txt) \
<(sort example_data/term2omim.py.txt)

*/