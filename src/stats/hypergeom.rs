//! Calculate the enrichment of a gene or disease in an `HpoSet` and the probability
//! of enrichment within the hypergeometric distribution.
//!
//! These methods are useful when you have clinical information of a patient
//! and want to see if the terms are enriched for some genes or diseases.
//!
//! # Examples
//!
//! ```
//! use hpo::term::InformationContentKind;
//! use hpo::{Ontology, HpoSet};
//! use hpo::term::HpoGroup;
//! use hpo::stats::hypergeom::gene_enrichment;
//!
//! fn clinial_info_set(ontology: &Ontology) -> HpoSet {
//! // ...
//! # let mut hpos = HpoGroup::new();
//! # hpos.insert(707u32);
//! # hpos.insert(12639u32);
//! # hpos.insert(12638u32);
//! # hpos.insert(818u32);
//! # hpos.insert(2715u32);
//! # HpoSet::new(ontology, hpos)
//! }
//!
//! let ontology = Ontology::from_binary("tests/example.hpo").unwrap();
//!
//! let patient_ci = clinial_info_set(&ontology);
//! let mut enrichments = gene_enrichment(&ontology, &patient_ci);
//!     
//! // the results are not sorted by default
//! enrichments.sort_unstable_by(|a, b| {
//!         a.pvalue().partial_cmp(&b.pvalue()).unwrap()
//!     });
//!
//! for gene in enrichments {
//!     println!("{}\t{}\t({})", gene.id(), gene.pvalue(), gene.enrichment());
//! }
//! ```
//! You can also use it to find genes that have a similar phenotype to a gene.
//! Here we are creating an `HpoSet` from all [`crate::HpoTerm`]s that are directly connected
//! to the gene `EZH2`. We're then checking which genes are enriched in the
//! directly and indirectly linked `HpoTerm`s.
//!
//! ```
//! use hpo::Ontology;
//! use hpo::{HpoSet, term::HpoGroup};
//! use hpo::stats::hypergeom::gene_enrichment;
//!
//! let ontology = Ontology::from_binary("tests/example.hpo").unwrap();
//!
//! let gene = ontology.gene_by_name("EZH2").unwrap();
//! let gene_hpo_set = gene.to_hpo_set(&ontology);
//!
//! let mut enrichments = gene_enrichment(&ontology, &gene_hpo_set);
//!
//! // the results are not sorted by default
//! enrichments.sort_by(|a, b| {
//!         a.pvalue().partial_cmp(&b.pvalue()).unwrap()
//! });
//! assert!(enrichments.first().unwrap().pvalue() < enrichments.last().unwrap().pvalue());
//! assert!(enrichments.first().unwrap().enrichment() > enrichments.last().unwrap().enrichment());
//! ```

mod disease;
mod gene;
pub use disease::disease_enrichment;
pub use gene::gene_enrichment;
