# HPO
This crate is currently a placeholder. I am trying to build a Rust version of [PyhPO](https://pypi.org/project/pyhpo/).

Since the library is not yet in a working state, this crate is empty.

If you find this project interesting, please get in touch, I could definitely need some help.

If you have another usecase for the `hpo` crate namespace and would like to use it, please let me know. I don't want to block the crate name if there are other or better use-cases for it.

## Basic requirements and use-cases
- Parse `obo` files and build an internal HPO ontology
- Parse annotation and metadata (genes and diseases)
- Connect HPO terms to genes and diseases
- Calculate the Information coefficient for all HPO terms
- Implement various similarity score algorithms
- Allow grouping of HPO terms into Sets, corresponding with the clinical information of a patient
- Implement similarity scores for sets
- Getting a single term from the ontology must be an O(1) operation
    - Access must be possible by ID and by name
- Traversing the hierarchy starting with one term must be fast
- Search functionality by term name (substring) (can be O(n) complex)
- Visualization of the ontology and relationshoips (e.g. via Graphviz)
- Implement Python bindings
- Be faster than PyHPO

### HPO-Term requirement
- Links to parent and child terms
- Links to genes and diseases

### Ontology requirements
- Allows iterating all terms
- Allows iterating all genes or diseases
- contains both HPO terms and annotations
- Provides methods to get terms and annotation in O(1) time

### Misc requirements
- API must provide an easy way to jump from one term to a related term. No need to pass around references to an arena and a term or something like it
- HPOTerms must be mutable initially and then can/should become immutable
    - We can only create the links to parents and children once all HPO terms are parsed

## Assumptions and Preconditions
- Ontology has around 20,000 terms (currently around 18.000)
- Around 20,000 genes and diseases (genes < 20,000, diseases ~10,000)
- Term ID can be easily converted to a unique integer
- Term ID is not continuous and there are huge gaps in between. The IDs are defined in batches
- The term ID is unique
- The term name is unique
- Diseases have a unique integer ID (non continuos)
- Genes have a unique integer ID (non continuos) and a unique name
- Focus only on human phenotypes, genes and diseases
- Currently disease sources are Omim, Orpha and Decipher

## Problems
- The ontology can have mixed types that are connected to each other
- HPO - HPO
- HPO - Gene
- HPO - Disease
- Gene - Disease

## Architecture and Design
### Data model of the Ontology
Setting up the ontology and individual terms brings some issues that complicate the overall setup. We must have a way to traverse the ontology from one term to another. This means that we must record all edges (connections) between the nodes (terms) in an efficient manner.
The `PyHPO` library utilzes Python's reference counting and aliasing functionality so that every term has references to all its direct children and direct parent terms. This allows efficient traversal, while at the same time allowing mutability of terms without sacrificing memory safety.
The mutability is important, because we can only link the references to parents and children after all terms are created. So we must rely on mutability during initialization. Once all terms and the ontology is fully built, we don't require mutability anymore.

Rust has a much stricter type system and it does not allow one object to have multiple references together with a mutable reference.
In [this very helpfully article](https://github.com/nrc/r4cppp/blob/master/graphs/README.md) two different options are suggesed. Using `Rc` (similar to Python's implementation) or using `UnsafeCell` and manage terms through an `arena`. Another option would similar to the arena approach: We could initialize all terms in a `Vec<HpoTerm>` and use the term-id (which is an integer) as the index of the record in the Vector. This would mean that the vector contains many empty records, but the Term-IDs are (almost) like an auto-increment number, so empty records would be sparse and should not be a problem for the memory.
However, this does not work if we want to work with pointers, since the memory location of the items changes when the vector is moved or resized.

Another option would be a dedicated graph data backend. The advantage would be that it includes several graph traversal and path finding algoriths that we could use. [petgraph](https://docs.rs/petgraph/latest/petgraph/index.html) looks very promising and could also print visual representations. In order to have constant time lookup, we would need a Hashmap that holds the HPO-Term-ID or HPO-Term-Name as key (`&str`) and the internal `NodeIndex` as value.

The main problem with a Graph is that we must have access to the Graph for all functions. So we could not call a method on HpoTerm or we would have to keep a reference to the Graph in every node (which sounds like it would cause a lot of issues with the borrow checker). Or every HPO term must keep references to all directly related HPO-terms.
```rust
struct HpoTerm {
    parents: Vec<&HpoTerm>,
    children: Vec<&HpoTerm>
}
```


I'm currently favoring the following idea:
Use a `typedArena::Arena` to store all `HpoTerm`s. Create a vector with `len() == max-id of HPO-terms` and fill with Option<&HpoTerm>

- Data model for Ontology
    - check existing ontology crates (e.g. `fastobo`), append functionality
    - check graph data structure implementations <>
    - provide some serialization of the data for faster startups
- obo file parser (input: File; output: Ontology)
- annotation metadata parser (input: File, ontology; output: Ontology)
- term similarity calculator (input: term-a, term-b; output: sim-score)
- HPO-set similarity calculator (input: term-a, term-b; output: sim-score)

Use <https://mermaid.live> for documentation


### API suggestions

#### Ontology
```rust
let ontology = some_function();

// iterate HPO terms
for term in ontology.hpos() {
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

#### HPO term
```rust
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