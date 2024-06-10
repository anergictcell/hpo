use rayon::prelude::*;
use std::io;
use std::io::Write;
use std::{env::Args, time::SystemTime};

use hpo::{
    annotations::{Disease, OmimDisease, OmimDiseaseId},
    similarity::{GraphIc, GroupSimilarity, StandardCombiner},
    term::HpoGroup,
    HpoSet, HpoTermId, Ontology,
};

/// Calculates the similarity score of two diseases
/// The two diseases are specified as OMIM-ID via CLI arguments
fn compare_two_diseases(
    ontology: &Ontology,
    sim: &GroupSimilarity<GraphIc, StandardCombiner>,
    mut args: Args,
) {
    let omimid_a = OmimDiseaseId::try_from(args.nth(1).unwrap().as_str())
        .expect("The first OmimID is invalid");
    let omimid_b = OmimDiseaseId::try_from(args.next().unwrap().as_str())
        .expect("The second OmimID is invalid");

    let omim_a = ontology
        .omim_disease(&omimid_a)
        .expect("The first OMIM Disease is not part of the Ontology");
    let omim_b = ontology
        .omim_disease(&omimid_b)
        .expect("The second OMIM Disease is not part of the Ontology");

    let set_a = omim_a.to_hpo_set(ontology);
    let set_b = omim_b.to_hpo_set(ontology);

    let res = sim.calculate(&set_a, &set_b);
    println!("Similarity is {res}");
}

/// Calculates the similarity score of a custom HPO-Set
/// to all OMIM diseases
/// The HPO-Set is specified as a comma separated list of HPO-Term-IDs
fn compare_custom_set_to_diseases(
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
    let mut results: Vec<(&OmimDisease, f32)> = ontology
        .omim_diseases()
        .par_bridge()
        .map(|disease| {
            let res = sim.calculate(&set_a, &disease.to_hpo_set(ontology));
            (disease, res)
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
fn cross_compare_diseases(
    ontology: &Ontology,
    sim: &GroupSimilarity<GraphIc, StandardCombiner>,
    num: usize,
) {
    let start = SystemTime::now();
    let results: Vec<(OmimDiseaseId, OmimDiseaseId, f32)> = ontology
        .omim_diseases()
        .par_bridge()
        .flat_map(|disease_a| {
            ontology
                .omim_diseases()
                .take(num)
                .map(|disease_b| {
                    let res = sim.calculate(
                        &disease_a.to_hpo_set(ontology),
                        &disease_b.to_hpo_set(ontology),
                    );
                    (*disease_a.id(), *disease_b.id(), res)
                })
                .collect::<Vec<(OmimDiseaseId, OmimDiseaseId, f32)>>()
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

fn compare_omim_to_orpha(
    ontology: &Ontology,
    sim: &GroupSimilarity<GraphIc, StandardCombiner>,
) {
    let omim: Vec<&OmimDisease> = ontology.omim_diseases().collect();

    let omim_names: Vec<&str> = omim.iter().map(|d| d.name()).collect();

    let omim_sets: Vec<HpoSet> = omim.iter().map(|d| d.to_hpo_set(ontology)).collect();

    println!("Orpha\\Omim\t{}", omim_names.join("\t"));

    ontology.orpha_diseases().take(100).par_bridge().for_each(|orpha| {
        let orpha_set = orpha.to_hpo_set(ontology);
        let mut row = orpha.name().to_string();
        for omim_set in omim_sets.iter() {
            row.push_str(&format!("\t{}", sim.calculate(&orpha_set, omim_set)));
        }
        println!("{row}");
    })
}

fn main() {
    let ontology = Ontology::from_binary("tests/ontology.hpo").unwrap();
    let combiner = StandardCombiner::FunSimAvg;
    let similarity = GraphIc::new(hpo::term::InformationContentKind::Omim);
    let sim = GroupSimilarity::new(combiner, similarity);

    let mut args = std::env::args();

    match args.len() {
        3 => compare_two_diseases(&ontology, &sim, args),
        2 => {
            let arg = args.nth(1).unwrap();
            if let Ok(num) = arg.parse::<usize>() {
                // integer provided, using disease x disease comparisons
                cross_compare_diseases(&ontology, &sim, num);
            } else if arg == "orpha" {
                let sim = GroupSimilarity::new(StandardCombiner::FunSimAvg, GraphIc::new(hpo::term::InformationContentKind::Gene));
                compare_omim_to_orpha(&ontology, &sim);
            }
            else {
                // List of HPO terms provided
                compare_custom_set_to_diseases(&ontology, &sim, arg);
            }
        }
        _ => {
            println!("Calculate similarities of OMIM diseases\n\n");
            println!("\
                There are 3 different options:\n\
                - Compare 2 diseases\n\
                    disease_similarity <OMIM ID> <OMIM ID>\n\
                    disease_similarity 618395 615368\n\n\
                - Compare all diseases to a custom HPO-Set\n\
                    disease_similarity <HPO-TERM-IDs>\n\
                    disease_similarity HP:0000750,HP:0000752,HP:0001249,HP:0007018,HP:0010818,HP:0011463\n\n\
                - Cross compare N diseases\n\
                    disease_similarity <NUMBER OF COMPARISONS>\n\
                    disease_similarity 20
                - Compare all OMIM to all ORPHA diseases\n\
                    disease_similarity orpha\n\
                    disease_similarity orpha
            ");
        }
    }
}
