# Changelog

All notable changes to this project will be documented in this file.

## Unreleased

## [0.11.1] - 2025-02-21

### Data

- Update to HPO 2025-01-16


## [0.11.0]

### Refactor

- [**breaking**] Add `Builder` to build Ontology. This removes several methods from `Ontology`.
- Remove aquamarine dependency, update docs


## [0.10.1]

### Refactor

- Add some missing methods for Orpha diseases


## [0.10.0]

### Feature

- Add Orphante diseases (`OrphaDisease`) to Ontology
- Filter gene and disease annotations in subontology based on association with phenotypes
- Add binary version 3
- Add new example ontology

### Documentation

- Change orders of methods in `Ontology` to clean up the documentation.

### Refactor

- Improve the OBO parser with better error handling
- [**breaking**] Add `Disease` trait that is needed to work with `OmimDisease` and `OrphaDisease`
- Update example ontology
- Update unit- and doctests to align with updated example ontology


## [0.9.1] - 2024-03-30

### Bugfix

- Fix the name of the `BMA` SimilarityCombiner.


## [0.9.0] - 2024-03-27

### Feature

- `Gene`s by default contain only direct `HpoTerm` associations, not transitive inherited ones.
- `Ontology::as_graohviz` method to generate graphviz data


## [0.8.3] - 2024-03-24

### Feature

- Add method to search for `OmimDisease`

### Bugfix

- Fix the JC similarity algorithm (see https://github.com/anergictcell/pyhpo/issues/20)

### Documentation

- Add a Changelog and a checklist for releases and patches


## [0.8.2] - 2024-03-09

### Data

- Update to HPO 2024-03-09

### Refactor

- Update dependencies


## [0.8.1] - 2023-06-25

### Feature

- Derive `Clone` for `Ontology`


## [0.8.0] - 2023-05-22

### Feature

- Add method to calculate hypergeometric enrichment of genes and diseases in HpoSets
- Add method to create dendogram clusters based on similarity

### Refactor

- Allow custom Similarity implementations to use Matrix


## [0.7.1] - 2023-04-27

### Refactor

- Derive `Debug` trait on more public structs


## [0.7.0] - 2023-04-22

### Feature

- New method to retrieve the shortest path between two HpoTerm
- Add modifier flag and categories of HpoTerm

### Refactor

- Use SmallVec for HpoGroup with default size 30
- Add more benchmarks
- Improve performance for adding, or-ing and comparing HpoGroups


## [0.6.3] - 2023-04-11

### Bugfix

- Fix issue parsing new HPO masterdata format


## [0.6.2] - 2023-04-05

### Bugfix

- Fix Subontology to not include all parents or children

### Refactor

- Add benchmark tests for Criterion


## [0.6.1] - 2023-03-30

### Documentation

- Add plenty of documentation


## [0.6.0] - 2023-03-18

### Feature

- Replace obsolete terms in an HpoSet
- allow different versions of binary masterdata

### Refactor

- add stricter clippy rules
- switch from `log` to `tracing`


## [0.5.0] - 2023-03-07

### Refactor

- clean up Similarity methods
- Simplify iterators across the full crate and add new ones


## [0.4.2] - 2023-02-11

### Feature

- new similarity method: Mutation


## [0.4.0] - 2023-02-04

### Feature

- Create a sub-ontology
- Calculate hypergeometric enrichment

### Bugfix

- Collecting into a HpoGroup will maintain order of the IDs internally
