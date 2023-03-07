use hpo::annotations::OmimDiseaseId;
use hpo::similarity::{
    Distance, GraphIc, InformationCoefficient, Jc, Lin, Relevance, Resnik, Similarity,
};
use hpo::term::InformationContentKind;
use hpo::Ontology;

/// We have to do this weird string reformatting, so that we can
/// compare the output to python. Python is printing the trailing 0,
/// e.g. `1.0` or `0.0`, but Rust isn't.
/// So we're adding the .0 for all int values.
fn round(n: f32, prec: usize) -> String {
    let factor = 10.0f32.powf(prec as f32);
    let res = (n * factor).round() / factor;
    if res % 1.0 == 0.0 {
        format!("{res}.0")
    } else {
        format!("{res}")
    }
}

fn scores_for_disease(omimid: OmimDiseaseId, ontology: &Ontology) {
    let resnik_o = Resnik::new(InformationContentKind::Omim);
    let resnik_g = Resnik::new(InformationContentKind::Gene);
    let lin_o = Lin::new(InformationContentKind::Omim);
    let lin_g = Lin::new(InformationContentKind::Gene);
    let jc_o = Jc::new(InformationContentKind::Omim);
    let jc_g = Jc::new(InformationContentKind::Gene);
    let rel_o = Relevance::new(InformationContentKind::Omim);
    let rel_g = Relevance::new(InformationContentKind::Gene);
    let ic_o = InformationCoefficient::new(InformationContentKind::Omim);
    let ic_g = InformationCoefficient::new(InformationContentKind::Gene);
    let graphic_o = GraphIc::new(InformationContentKind::Omim);
    let graphic_g = GraphIc::new(InformationContentKind::Gene);
    let dist_o = Distance::new();

    let disease = ontology.omim_disease(&omimid).unwrap();
    let terms = disease.to_hpo_set(ontology);
    for term1 in &terms {
        for term2 in &terms {
            println!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                term1.id(),
                term2.id(),
                round(resnik_o.calculate(&term1, &term2), 4),
                round(resnik_g.calculate(&term1, &term2), 4),
                round(lin_o.calculate(&term1, &term2), 4),
                round(lin_g.calculate(&term1, &term2), 4),
                round(jc_o.calculate(&term1, &term2), 4),
                round(jc_g.calculate(&term1, &term2), 4),
                round(rel_o.calculate(&term1, &term2), 4),
                round(rel_g.calculate(&term1, &term2), 4),
                round(ic_o.calculate(&term1, &term2), 4),
                round(ic_g.calculate(&term1, &term2), 4),
                round(graphic_o.calculate(&term1, &term2), 4),
                round(graphic_g.calculate(&term1, &term2), 4),
                round(dist_o.calculate(&term1, &term2), 4),
            );
        }
    }
}

fn main() {
    let ontology = Ontology::from_binary("tests/ontology.hpo").unwrap();

    scores_for_disease("300486".try_into().unwrap(), &ontology);
}

/*
from pyhpo import Ontology
from pyhpo.annotations import Omim

_ = Ontology()

disease = Omim.get(300486)
# OmimDisease(id=300486, name='Mental retardation, X-linked, with cerebellar hypoplasia and distinctive facial appearance', hpo={256, 3593, 1290, 1419, 400, 2066, 276, 7065, 28, 6951, 1320, 1321, 46, 303, 30260, 54, 1344, 31936, 448, 322, 2119, 336, 11220, 2007, 601, 219, 733, 1249, 1250, 486, 742, 2280, 744, 490, 1257, 1263, 752, 3189, 25336, 639}, negative_hpo=set(), diseasetype='Omim')

terms = [Ontology[idx] for idx in disease.hpo]

with open("similarity_batch.py.txt", "w") as fh:
    for term1 in terms:
        for term2 in terms:
            ic = term1.similarity_score(term2, kind="omim", method="graphic")
            _ = fh.write((
                    "{}\t{}\t{}\n".format(
                        term1.id,
                        term2.id,
                        "\t".join([
                            "{}".format(round(term1.similarity_score(term2, kind="omim", method="resnik"), 4)),
                            "{}".format(round(term1.similarity_score(term2, kind="gene", method="resnik"), 4)),
                            "{}".format(round(term1.similarity_score(term2, kind="omim", method="lin"), 4)),
                            "{}".format(round(term1.similarity_score(term2, kind="gene", method="lin"), 4)),
                            "{}".format(round(term1.similarity_score(term2, kind="omim", method="jc2"), 4)),
                            "{}".format(round(term1.similarity_score(term2, kind="gene", method="jc2"), 4)),
                            "{}".format(round(term1.similarity_score(term2, kind="omim", method="rel"), 4)),
                            "{}".format(round(term1.similarity_score(term2, kind="gene", method="rel"), 4)),
                            "{}".format(round(term1.similarity_score(term2, kind="omim", method="ic"), 4)),
                            "{}".format(round(term1.similarity_score(term2, kind="gene", method="ic"), 4)),
                            "{}".format(round(term1.similarity_score(term2, kind="omim", method="graphic"), 4)),
                            "{}".format(round(term1.similarity_score(term2, kind="gene", method="graphic"), 4)),
                            "{}".format(round(term1.similarity_score(term2, kind="omim", method="dist"), 4)),
                        ])
                    )
                ))
*/
