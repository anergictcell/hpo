use hpo::{annotations::OmimDiseaseId, term::HpoGroup, HpoSet, Ontology};

use std::time::SystemTime;

fn main() {
    // simple_logger::init_with_level(log::Level::Debug).unwrap();
    let ontology = Ontology::from_binary("tests/ontology.hpo").unwrap();

    let mut args = std::env::args();
    if args.len() < 2 {
        panic!("You must provide at least 1 argument")
    }
    let hpo_set = match args.nth(1).unwrap().as_str() {
        "gene" => ontology
            .gene_by_name(&args.next().unwrap())
            .expect("Invalid gene symbol")
            .to_hpo_set(&ontology),
        "disease" => ontology
            .omim_disease(&OmimDiseaseId::try_from(args.next().unwrap().as_str()).unwrap())
            .expect("Invalid Omim disease")
            .to_hpo_set(&ontology),
        terms => {
            let mut hpos = HpoGroup::new();
            for id in terms.split(',') {
                hpos.insert(id.parse::<u32>().unwrap().into());
            }
            HpoSet::new(&ontology, hpos)
        }
    };

    let start = SystemTime::now();
    let mut gene_enrichments = hpo::stats::hypergeom::gene_enrichment(&ontology, &hpo_set);
    let end = SystemTime::now();
    let duration = end.duration_since(start).unwrap();

    gene_enrichments.sort_by(|a, b| {
        a.pvalue()
            .partial_cmp(&b.pvalue())
            .expect("nan must not appear as enrichment")
    });
    println!("### GENES ###");
    for gene in &gene_enrichments[0..10] {
        println!(
            "{}\t{:e}\t({})",
            ontology.gene(gene.id()).unwrap().name(),
            gene.pvalue(),
            gene.enrichment()
        );
    }
    println!(
        "It took {} milliseconds for {} genes",
        duration.as_millis(),
        gene_enrichments.len()
    );

    let start = SystemTime::now();
    let mut disease_enrichments = hpo::stats::hypergeom::disease_enrichment(&ontology, &hpo_set);
    let end = SystemTime::now();
    let duration = end.duration_since(start).unwrap();

    disease_enrichments.sort_by(|a, b| {
        a.pvalue()
            .partial_cmp(&b.pvalue())
            .expect("nan must not appear as enrichment")
    });
    println!("\n\n### DISEASES ###");
    for disease in &disease_enrichments[0..10] {
        println!(
            "{}\t{:e}\t({})",
            ontology.omim_disease(disease.id()).unwrap().name(),
            disease.pvalue(),
            disease.enrichment()
        );
    }
    println!(
        "It took {} milliseconds for {} diseases",
        duration.as_millis(),
        disease_enrichments.len()
    );

    println!(
        "\nTerms: {}\nTotal gene: {}\nTotal diseases: {}",
        hpo_set.len(),
        ontology.genes().count(),
        ontology.omim_diseases().count()
    );

    println!("\n\n## Matches");
    println!("### Diseases");

    if let Some(filter) = args.next() {
        for d in disease_enrichments.iter().filter_map(|d| {
            let dn = ontology.omim_disease(d.id()).unwrap().name();
            if dn.to_lowercase().contains(&filter.to_lowercase()) {
                Some(format!("{dn}\t{}\t{}", d.pvalue(), d.enrichment()))
            } else {
                None
            }
        }) {
            println!("{}", d)
        }

        println!("### Genes");
        for d in gene_enrichments.iter().filter_map(|d| {
            let dn = ontology.gene(d.id()).unwrap().name();
            if dn.to_lowercase().contains(&filter.to_lowercase()) {
                Some(format!("{dn}\t{}\t{}", d.pvalue(), d.enrichment()))
            } else {
                None
            }
        }) {
            println!("{}", d)
        }
    }
}
