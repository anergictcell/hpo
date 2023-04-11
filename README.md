# HPO

This library is a Rust implementation of [PyHPO](https://pypi.org/project/pyhpo/).

## What is this?

HPO, the [Human Phenotype Ontology](https://hpo.jax.org/app/) is a standard vocabulary of phenotypic abnormalities in human diseases. It is an Ontology, so all terms are connected to each other, similar to a directed graph.  

This library provides convenient APIs to work with the ontology. The main goals are to compare terms - or sets of terms - to each other and run statistics for enrichment analysis.

### Features
- Calculate the similarity of HPO terms
- Calculate the similarity of multiple sets of HPO terms (e.g. a patient's clinical information)
- Enrichment analysis of genes and diseases in sets of HPO terms
- Compare different HPO versions
- Graph based analysis of the ontology
- Completely written in Rust, so it's **🚀blazingly fast🚀**<sup>TM</sup> ([Benchmarks](#benchmarks))

## What is the current state?

The library is pretty much feature-complete, at least for my use-cases. If you have any feature-requests, please open an Issue or get in touch. I'm very much interested in getting feedback and new ideas what to improve.

The API is mostly stable, but I might refactor some parts a bit for easier use and performance gain.

If you find this project interesting and want to contribute, please get in touch, I could definitely need some help.

## Documentation
The public API is fully documented on [`docs.rs`](https://docs.rs/hpo/latest/hpo/)

The main structs used in `hpo` are:
- The [`Ontology`](https://docs.rs/hpo/latest/hpo/struct.Ontology.html) is the main struct and entrypoint in `hpo`.
- [`HpoTerm`](https://docs.rs/hpo/latest/hpo/term/struct.HpoTerm.html) represents a single HPO term and contains plenty of functionality around them.
- [`HpoSet`](https://docs.rs/hpo/latest/hpo/struct.HpoSet.html) is a collection of `HpoTerm`s, like a patient's clinical information.
- [`Gene`](https://docs.rs/hpo/latest/hpo/annotations/struct.Gene.html) represents a single gene, including information about associated `HpoTerm`s.
- [`OmimDisease`](https://docs.rs/hpo/latest/hpo/annotations/struct.OmimDisease.html) represents a single OMIM-diseases, including information about associated `HpoTerm`s.

The most relevant modules are:
- [`annotations`](https://docs.rs/hpo/latest/hpo/annotations/index.html) contains the `Gene` and `OmimDisease` structs, and some related important types.
- [`similarity`](https://docs.rs/hpo/latest/hpo/similarity/index.html) contains structs and helper functions for similarity comparisons for `HpoTerm` and `HpoSet`.
- [`stats`](https://docs.rs/hpo/latest/hpo/stats/index.html) contains functions to calculate the hypergeometric enrichment score of genes or diseases.


## Examples
Some (more or less random) examples are included in the [`examples` folder](https://github.com/anergictcell/hpo/tree/main/examples).

### Ontology
```rust
use hpo::{Ontology, HpoTermId};
use hpo::annotations::{GeneId, OmimDiseaseId};

fn example() {
    let ontology = Ontology::from_standard("/path/to/master-data/").unwrap();

    // iterate HPO terms
    for term in &ontology {
        // do something with term
    }

    // iterate Genes
    for gene in ontology.genes() {
        // do something with gene
    }

    // iterate omim diseases
    for disease in ontology.omim_diseases() {
        // do something with disease
    }

    // get a single HPO term using HPO ID
    let hpo_id = HpoTermId::try_from("HP:0000123").unwrap();
    let term = ontology.hpo(hpo_id);

    // get a single HPO term using `u32` part of HPO ID
    let hpo_id = HpoTermId::from(123u32);
    let term = ontology.hpo(hpo_id);

    // get a single Omim disease
    let disease_id = OmimDiseaseId::from(12345u32);
    let disease = ontology.omim_disease(&disease_id);

    // get a single Gene
    let hgnc_id = GeneId::from(12345u32);
    let gene = ontology.gene(&hgnc_id);

    // get a single Gene by its symbol
    let gene = ontology.gene_by_name("GBA");

}
```

### HPO term
```rust
use hpo::Ontology;

fn example() {
    let ontology = Ontology::from_binary("/path/to/binary.hpo").unwrap();

    let term = ontology.hpo(123u32.into()).unwrap();

    assert_eq!("Abnormality of the nervous system", term.name());
    assert_eq!("HP:000123".to_string(), term.id().to_string());

    // iterate all parents
    for p in term.parents() {
        println!("{}", p.name())
    }

    // iterate all children
    for p in term.children() {
        println!("{}", p.name())
    }

    let term2 = ontology.hpo(1u32.into()).unwrap();

    assert!(term2.parent_of(&term));
    assert!(term.child_of(&term2));
}
```

### Similarity
```rust
use hpo::Ontology;
use hpo::similarity::GraphIc;
use hpo::term::InformationContentKind;

fn example() {
    let ontology = Ontology::from_binary("/path/to/binary.hpo").unwrap();
    let term1 = ontology.hpo(123u32.into()).unwrap();
    let term2 = ontology.hpo(1u32.into()).unwrap();

    let ic = GraphIc::new(InformationContentKind::Omim);
    let similarity = term1.similarity_score(&term2, &ic);
}
```

### Enrichment
Identify which genes (or diseases) are enriched in a set of HpoTerms, e.g. in
the clinical information of a patient or patient cohort

```rust
use hpo::Ontology;
use hpo::{HpoSet, term::HpoGroup};
use hpo::stats::hypergeom::gene_enrichment;

fn example() {
    let ontology = Ontology::from_binary("/path/to/binary.hpo").unwrap();

    let mut hpos = HpoGroup::new();
    hpos.insert(2943u32.into());
    hpos.insert(8458u32.into());
    hpos.insert(100884u32.into());
    hpos.insert(2944u32.into());
    hpos.insert(2751u32.into());
    let patient_ci = HpoSet::new(&ontology, hpos);

    let mut enrichments = gene_enrichment(&ontology, &patient_ci);

    // the results are not sorted by default
    enrichments.sort_by(|a, b| {
        a.pvalue().partial_cmp(&b.pvalue()).unwrap()
    });

    for gene in enrichments {
        println!("{}\t{}\t({})", gene.id(), gene.pvalue(), gene.enrichment());
    }
}
```

## Benchmarks
As the saying goes: "Make it work, make it good, make it fast". The *work* and *good* parts are realized in [PyHPO](https://pypi.org/project/pyhpo/). And even though I tried my best to make it *fast*, I was still hungry for more. So I started developing the `hpo` Rust library in December 2022. Even without micro-benchmarking and tuning performance as much as I did for `PyHPO`, `hpo` is indeed much much faster already now.

The below benchmarks were run non scientificially and your mileage may vary. I used a MacBook Air M1, `rustc 1.68.0`, `Python 3.9` and `/usr/bin/time` for timing.

| Benchmark | `PyHPO` | `hpo` (single-threaded) | `hpo` (multi-threaded) |
| --------- | ----- | --- | --- |
| Read and Parse Ontology | 6.4 s | 0.3 s | 0.3 s |
| Similarity of 17,245 x 1,000 terms | 98.5 s | 6.2 s | 1.7 s |
| Similarity of GBA1 to all Diseases | 380 s | 27.3 s | 6.1 s |
| Disease enrichment in all Genes | 11.8 s | 0.6 s | 0.3 s |
| Common ancestors of 17,245 x 10,000 terms | 200.6 s | 12.1 | 2.8 |



## Technical design
There is some info about the plans for the implementation in the [Technical Design document](https://github.com/anergictcell/hpo/blob/main/TechnicalDesign.md)
