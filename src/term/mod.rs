//! [`HpoTerm`]s are the main building block of the Ontology.
//!  Each term is a descendent (child) of at least one other term
//! (except for the root term `HP:0000001 | All`). The relationship
//! is modeled bi-drectionally in `hpo`, so that every term also has
//! one or several `children` (except for all leaf terms).
//!
//! At the same time, terms are associated with diseases by describing a part
//! of their phenotypical representation. They are also associated with genes
//! when the gene is relevant for the phenotype or associated disease. The
//! association to genes and diseases are inversely "inherited" to a term's parent
//! terms (and subsequently their parents as well).
//!
//! Each term is identified by a unique [`HpoTermId`].
//!
use crate::Ontology;
use core::fmt::Debug;

pub mod group;
mod hpoterm;
mod hpotermid;
mod information_content;
pub(crate) mod internal;

pub use group::HpoGroup;
pub use hpoterm::HpoTerm;
pub use hpotermid::HpoTermId;
pub use information_content::{InformationContent, InformationContentKind};

/// [`HpoTerm`] iterator
pub struct Iter<'a> {
    group_iter: group::Iter<'a>,
    ontology: &'a Ontology,
}

impl<'a> Iter<'a> {
    pub(crate) fn new(group_iter: group::Iter<'a>, ontology: &'a Ontology) -> Self {
        Self {
            group_iter,
            ontology,
        }
    }
}

impl<'a> Iterator for Iter<'a> {
    type Item = HpoTerm<'a>;
    fn next(&mut self) -> Option<Self::Item> {
        match self.group_iter.next() {
            Some(term) => {
                let term = self
                    .ontology
                    .get(term)
                    .unwrap_or_else(|| panic!("Invalid HPO-Term: {term}"));
                Some(HpoTerm::new(self.ontology, term))
            }
            None => None,
        }
    }
}

impl Debug for Iter<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Iter<HpoTerm>")
    }
}
