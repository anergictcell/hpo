//! Contains structs related to HPO-Terms

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
    group: HpoTermIds<'a>,
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
