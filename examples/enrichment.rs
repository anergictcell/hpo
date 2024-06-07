use std::{env::Args, process};

use hpo::{
    annotations::{Disease, OmimDiseaseId},
    term::HpoGroup,
    HpoResult, HpoSet, HpoTermId, Ontology,
};

/// Tries to parse an HpoTermId from a string `HP:0007768` or a `u32`
fn id_from_freetext(value: &str) -> HpoResult<HpoTermId> {
    if let Ok(i) = value.parse::<u32>() {
        Ok(HpoTermId::from(i))
    } else {
        HpoTermId::try_from(value)
    }
}

/// Prints genes that are enriched in the hposet
fn gene_enrichments(ontology: &Ontology, hposet: &HpoSet, output_len: usize) {
    let mut gene_enrichments = hpo::stats::hypergeom::gene_enrichment(ontology, hposet);

    gene_enrichments.sort_by(|a, b| {
        a.pvalue()
            .partial_cmp(&b.pvalue())
            .expect("nan must not appear as enrichment")
    });
    println!("### GENES ###");
    for gene in &gene_enrichments[0..std::cmp::min(output_len, gene_enrichments.len())] {
        println!(
            "{}\t{:e}\t({})",
            ontology.gene(gene.id()).unwrap().name(),
            gene.pvalue(),
            gene.enrichment()
        );
    }
}

/// Prints diseases that are enriched in the hposet
fn disease_enrichments(ontology: &Ontology, hposet: &HpoSet, output_len: usize) {
    let mut disease_enrichments = hpo::stats::hypergeom::disease_enrichment(ontology, hposet);

    disease_enrichments.sort_by(|a, b| {
        a.pvalue()
            .partial_cmp(&b.pvalue())
            .expect("nan must not appear as enrichment")
    });
    println!("\n\n### DISEASES ###");
    for disease in &disease_enrichments[0..std::cmp::min(output_len, disease_enrichments.len())] {
        println!(
            "{}\t{:e}\t({})",
            ontology.omim_disease(disease.id()).unwrap().name(),
            disease.pvalue(),
            disease.enrichment()
        );
    }
}

/// Parses an HpoSet from a gene symbol, disease ID or a list of HpoTermIds
fn hposet<'a>(args: &mut Args, ontology: &'a Ontology) -> HpoSet<'a> {
    match args.nth(1).unwrap().as_str() {
        "gene" => ontology
            .gene_by_name(&args.next().unwrap())
            .expect("Invalid gene symbol")
            .to_hpo_set(ontology),
        "disease" => ontology
            .omim_disease(&OmimDiseaseId::try_from(args.next().unwrap().as_str()).unwrap())
            .expect("Invalid Omim disease")
            .to_hpo_set(ontology),
        terms => {
            let mut hpos = HpoGroup::new();
            for id in terms.split(',') {
                hpos.insert(id_from_freetext(id).expect("invalid HPO ID"));
            }
            HpoSet::new(ontology, hpos)
        }
    }
}

fn main() {
    let ontology = Ontology::from_binary("tests/ontology.hpo").unwrap();

    let mut args = std::env::args();
    if args.len() < 2 {
        println!("Show enriched genes and diseases\n\n");
        println!("Usage\nenrichment gene GBA1 <N RESULTS>");
        println!("\nenrichment disease 607417 <N RESULTS>");
        println!("\nenrichment HP:0007768,HP:0000631 <N RESULTS>\n");
        process::exit(1)
    }

    let hpo_set = hposet(&mut args, &ontology);

    let output_len = args
        .next()
        .map(|arg| arg.parse::<usize>().unwrap_or(10))
        .unwrap_or(10);

    gene_enrichments(&ontology, &hpo_set, output_len);
    disease_enrichments(&ontology, &hpo_set, output_len);

    println!(
        "\nTerms: {}\nTotal gene: {}\nTotal diseases: {}",
        hpo_set.len(),
        ontology.genes().count(),
        ontology.omim_diseases().count()
    );
}
