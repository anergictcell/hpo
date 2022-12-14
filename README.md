# HPO

This library is a Rust implementation of [PyHPO](https://pypi.org/project/pyhpo/).

## What is this?

HPO, the [Human Phenotype Ontology](https://hpo.jax.org/app/) is a standard vocabulary of phenotypic abnormalities in human diseases. It is an Ontology, so all terms are connected to each other, similar to a directed graph.  
This library provides convenient APIs to work with the ontology. The main goals are to compare terms to each other and compare group of terms to each other.

For example, the terms `"Migraine without aura"` and `"Migraine with aura"` are more similar to each other than `"Migraine"` and `"Seizure"`. To add more complexity, patients usually have more than one phenotypical abnormality. So in order to compare two patients to each other, we must cross-compare all individual terms. Eventually we might want to cluster hundreds or thousands of patients based on phenotypical similarity to predict diseases based on phenotypes and run statistical analyses.

The [PyHPO](https://pypi.org/project/pyhpo/) Python library provides functionality for these comparisons, providing several different similarity and grouping algorithms. However, since its written in Python it is very slow. Unfortunately the design of PyHPO does not allow multithreading or parallel processing, which makes scaling rather difficult.

I want to overcome these limitations here with a Rust library.

There is some info about the plans for the implementation in the [Technical Design document](https://github.com/anergictcell/hpo/blob/main/TechnicalDesign.md)


If you find this project interesting and want to contribute, please get in touch, I could definitely need some help. The code is not yet well documented and does not yet have many tests. At the moment, I'm primarily trying to get a working PoC and use the `examples` to test and compare outputs to `PyHPO`.

## What is the current state?

At the moment, this library provides most of the functionality of `PyHPO` and it does so much much faster (**Blazingly fast**). For example, to calculate the `GraphIC` similarity for 400 x 400 terms, `PyHPO` runs for about 40 seconds on my MacBook Air M1. This Rust based `hpo` library finishes in less than 1 second. I can run a pairwise comparison of all 17,059 terms to each other in 13 seconds. `PyHPO` would need several hours for the same task.
`hpo` also allows multithreading, e.g. using rayon.

You can check out some examples, including benchmarks, in the `examples` folder. I sometimes also include the corresponding Python code. As with all benchmarks, your mileage may vary, depending on your computer. But the overall trend stays the same.

## API suggestions

At the moment, not all functionality in here is working, but most of it is. Check out the `examples` folder in the Github repository for some ways to use it.

### Ontology
```rust
use hpo::Ontology;
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

    // get a single HPO term
    let term = ontology.hpo("HP:0000123".try_into().unwrap());

    // get a single Gene
    let hgnc_id = GeneId::try_from("12345").unwrap();
    let gene = ontology.gene(&hgnc_id);

    // get a single Omim disease
    let disease_id = OmimDiseaseId::try_from("12345").unwrap();
    let disease = ontology.omim_disease(&disease_id);
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