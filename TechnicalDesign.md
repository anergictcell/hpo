# Technical Software Description

## Basic requirements and use-cases

- Parse `obo` files and build an internal HPO ontology
- Parse annotation and metadata (genes and diseases)
- Connect `HPO terms` to genes and diseases
- Calculate the Information content (IC) for all `HPO terms`
- Implement various similarity score algorithms
- Allow grouping of `HPO terms` into `Sets`, corresponding with the clinical information of a patient
- Implement similarity scores for `Sets`
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
- is the term obsolete or even replaces
- List of all direct and indirect parents

### Ontology requirements

- Allows iterating all `HPO terms`
- Allows iterating all genes or diseases
- contains both `HPO terms` and annotations
- Provides methods to get `HPO terms` and annotation in O(1) time

### Misc requirements

- API must provide an easy way to jump from one term to a related term
- HPOTerms must be mutable initially and then become immutable
    - We can only create the links to parents and children once all HPO terms are parsed

## Assumptions and Preconditions

- Ontology has around 20,000 terms (currently around 18.000).
- Around 20,000 genes and diseases (genes < 20,000, diseases ~10,000).
- Term ID can be easily converted to a unique integer.
- Term ID is not continuous and there are huge gaps in between.
- The term ID is unique for every term.
- The term name is unique.
- Diseases have a unique integer ID (non continuos).
- Genes have a unique integer ID (non continuos) and a unique name.
- Focus only on human phenotypes, genes and diseases.
- Currently disease sources are Omim, Orpha and Decipher.
- Some terms are obsolete and only kept for backwards compatibility. They are not connected to others. Some obsolete terms indicate a replacement term.
- The ontology can have mixed types that are connected to each other
    - HPO - HPO
    - HPO - Gene
    - HPO - Disease
    - Gene - Disease

## Architecture and Design

### Data model of the Ontology
The structure of the ontology and the connections between the different entities is quite complex and I could not find a good matching data model in the standard library or as a 3rd party crate.
The most critical functionality is the traversal of the ontology from one term to another and identify common ancestor terms. This means that we must record all edges (connections) between the nodes (terms) in an efficient manner. A term should look a bit like this:

```rust
struct HpoTerm {
    parents: Vec<&HpoTerm>,
    children: Vec<&HpoTerm>,
    // ... other fields
}
```

The `PyHPO` library utilzes Python's reference counting and aliasing functionality so that every term has references to all its direct children and parent terms. This allows efficient traversal, while at the same time allowing mutability of terms without sacrificing memory safety.
The mutability is important, because we must modify HPO terms during the creation of the ontology. For the creation, we first must initialize all terms. Only once all terms are present, can we link them to each other. So we must rely on mutability during initialization. Once all terms and the ontology is fully built, we don't require mutability anymore.

Rust has a much stricter type system and it does not allow one object to have multiple references together with a mutable reference.
In [this very helpfully article](https://github.com/nrc/r4cppp/blob/master/graphs/README.md) two different options are suggesed. Using `Rc` (similar to Python's implementation) or using `UnsafeCell` and manage terms through an `arena`.  
Another option would be a dedicated graph data backend. The advantage would be that it includes several graph traversal and path finding algoriths that we could use. [petgraph](https://docs.rs/petgraph/latest/petgraph/index.html) looks very promising and could also print visual representations. In order to have constant time lookup, we would need a Hashmap that holds the HPO-Term-ID or HPO-Term-Name as key (`&str`) and the internal `NodeIndex` as value.


### Current approach

I am currently building a PoC with the `arena`-like approach. All `HPO terms` are stored in a `Vec<HpoTerm` and we use a lookup table to jump from the `HpoTermId` (derived from the HPO-ID) to the actual `HpoTerm` in the `Vec`. Initially I tried using a `HashMap` as lookup table, but found out that a large `Vec` is much faster.

This means the Ontology consists of one (small) `Vec` which contains every HPO terms in no particular order (let's call it `Term-Vec`) and one large vector which contains one record for every possible HPO-Term ID (lets call it `ID-Vec`). Since every HPO-Term ID can be converted into an integer in the range of 1 - 10,000,000 this vector has a length of 10 Mio. Every record is a usize which indicates the index of the `Term-Vec`. Due to this we can very efficiently retrieve an HPO-Term from its ID:

1. Convert ID (`HP:000123`) to its integer `123` representation
2. get the record at index `123` from the `ID-Vec` (e.g: 17)
3. Use that record as index for the `Term-Vec` to retrieve the actual HPO term