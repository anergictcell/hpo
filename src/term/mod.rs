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
//! 'hpo' supports many operations that work on sets of `HpoTerms`, they differ
//! in their implementation details. The two most common are:
//! - [`HpoTerms`]: An iterator of [`HpoTerm`] structs. In most cases this
//!   represents a set, i.e. every term is unique, although that is not guaranteed
//! - [`HpoGroup`]: A set of [`HpoTermId`] implements various functionality
//! - [`HpoTermIds`]: An iterator of [`HpoTermId`], which can be converted from/into an [`HpoGroup`]
//!
use crate::Ontology;
use core::fmt::Debug;

mod group;
mod hpoterm;
mod hpotermid;
mod information_content;
pub(crate) mod internal;

pub use group::{HpoGroup, HpoTermIds};
pub use hpoterm::HpoTerm;
pub use hpotermid::HpoTermId;
pub use information_content::{InformationContent, InformationContentKind};

/// A set of parent [`HpoTermId`]s
///
/// It uses [`HpoGroup`] underneath
pub type HpoParents = HpoGroup;

/// A set of child [`HpoTermId`]s
///
/// It uses [`HpoGroup`] underneath
pub type HpoChildren = HpoGroup;

/// Iterate [`HpoTerm`]s
///
/// This struct creates [`HpoTerm`]s from a reference to [`HpoGroup`]
pub struct HpoTerms<'a> {
    ontology: &'a Ontology,
    group: HpoTermIds<std::slice::Iter<'a, HpoTermId>>,
}

impl<'a> HpoTerms<'a> {
    /// Returns a new [`HpoTerms`]
    #[must_use]
    pub fn new(group: &'a HpoGroup, ontology: &'a Ontology) -> Self {
        HpoTerms {
            group: group.iter(),
            ontology,
        }
    }

    pub fn ontology(&self) -> &Ontology {
        self.ontology
    }
}

impl<'a> Iterator for HpoTerms<'a> {
    type Item = HpoTerm<'a>;
    fn next(&mut self) -> Option<Self::Item> {
        match self.group.next() {
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

impl Debug for HpoTerms<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "HpoTermIterator")
    }
}
