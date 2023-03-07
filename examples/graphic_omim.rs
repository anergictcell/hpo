use hpo::annotations::OmimDiseaseId;
use hpo::similarity::Builtins;
use hpo::term::InformationContentKind;
use hpo::Ontology;

fn scores_for_disease(omimid: OmimDiseaseId, ontology: &Ontology) {
    let ic = Builtins::GraphIc(InformationContentKind::Omim);

    let disease = ontology.omim_disease(&omimid).unwrap();
    println!("Using Disease [{}] {}", disease.id(), disease.name());
    let terms = disease.to_hpo_set(ontology);
    for term1 in &terms {
        for term2 in &terms {
            let overlap = term1.similarity_score(&term2, &ic);
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
    let ontology = Ontology::from_standard("./example_data/").unwrap();

    scores_for_disease("300486".try_into().unwrap(), &ontology);
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
