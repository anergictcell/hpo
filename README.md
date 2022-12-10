# HPO

This crate is a draft for a Rust implementation of [PyHPO](https://pypi.org/project/pyhpo/).

> :warning: **Warning:** The library is a work in progress and I don't recommend using it at this stage.

If you find this project interesting and want to contribute, please get in touch, I could definitely need some help. The code is not yet well documented and does not yet have many tests. At the moment, I'm primarily trying to get a working PoC. Once I'm there, I will adjust many method names and functionality and add more documentation and tests. The library does not contain any error handling and uses `unwrap` a lot - I plan to change this once I am ready to stabilize the overall API a bit more.

If you have another usecase for the `hpo` crate namespace and would like to use it, please let me know. I don't want to block the crate name if there are better use-cases for it.

## What is this?

HPO, the [Human Phenotype Ontology](https://hpo.jax.org/app/) is a standard vocabulary of phenotypic abnormalities in human diseases. It is an Ontology, so all terms are connected to each other, similar to a directed graph.  
This crate should give some convinient APIs to work with the ontology. The main goals are to compare terms to each other and also compare group of terms to each other.
For example, the terms `Migraine without aura` and `Migraine with aura` are more similar to each other than `Migraine` and `Seizure`. To add more complexity, patients usually have more than one phenotypical abnormality. So in order to compare two patients to each other, we must cross-compare all individual terms. Eventually we might want to cluster hundreds or thousands of patients based on phenotypical similarity.

The [PyHPO](https://pypi.org/project/pyhpo/) Python library provides functionality for these comparisons, providing several different similarity and grouping algorithms. However, since its written in Python it is rather slow. Unfortunately the design of PyHPO does not allow multithreading or parallel processing, which makes scaling rather difficult.

I want to overcome these limitations here with a Rust library.

There is some info about the plans for the implementation in the [Technical Design document](https://github.com/anergictcell/hpo/blob/main/TechnicalDesign.md)


## API suggestions
At the moment, not all functionality in here is working, but most of it is. Check out the `examples` folder in the Github repository for some ways to use it.

### Ontology
```ignore
let ontology = some_function();

// iterate HPO terms
for term in ontology {
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
let term = ontology.hpo("HP:0000123");
let term = ontology.hpo_by_name("Abnormality of the nervous system");

// search HPO terms
for res in ontology.search_hpo("Abnormality") {
    println!("{}", res.name());
}

// get a single Gene
let hgnc_id = 12345;
let gene = ontology.gene(hgnc_id);
let gene = ontology.gene_by_name("EZH2");

// get a single Omim disease
let disease_id = 12345;
let disease = ontology.omim_disease(disease_id);
let disease = ontology.omim_by_name("Gaucher");
```

### HPO term
```ignore
let term = some_function();

assert_eq!("Abnormality of the nervous system", term.name());
assert_eq!("HP:000123", term.id());

// iterate all parents
for p in term.parents() {
    println!("{}", p.name())
}

// iterate all children
for p in term.children() {
    println!("{}", p.name())
}

let term2 = some_other_function()

assert!(term.parent_of(term2));
assert!(term2.child_of(term));

assert!(!term2.parent_of(term));
assert!(!term.child_of(term2));
```

### Similarity
```ignore
let term = some_function();
let term2 = some_other_function()
let ic = GraphIc::new(hpo::InformationContentKind::Omim);
let similarity = term1.similarity_score(&term2, &ic);
```