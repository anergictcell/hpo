use core::fmt::Debug;
use crate::Ontology;

mod group;
mod hpoterm;
mod hpotermid;
mod information_content;
pub (crate) mod internal;
mod parents;

pub use group::{HpoGroup, HpoGroupIterator};
pub use hpoterm::HpoTerm;
pub use hpotermid::HpoTermId;
pub use information_content::{InformationContent, InformationContentKind};

pub type HpoParents = HpoGroup;
pub type HpoChildren = HpoGroup;


pub struct HpoTermIterator<'a> {
    ontology: &'a Ontology,
    parents: HpoGroupIterator<'a>,
}

impl<'a> HpoTermIterator<'a> {
    pub fn new(parents: &'a HpoGroup, ontology: &'a Ontology) -> Self {
        HpoTermIterator {
            parents: parents.iter(),
            ontology,
        }
    }
}

impl<'a> std::iter::Iterator for HpoTermIterator<'a> {
    type Item = HpoTerm<'a>;
    fn next(&mut self) -> Option<Self::Item> {
        match self.parents.next() {
            Some(term) => {
                let term = self.ontology.get(term).unwrap_or_else(|| panic!("Invalid HPO-Term: {}", term));
                Some(HpoTerm::new(self.ontology, term))
            }
            None => None,
        }
    }
}

impl Debug for HpoTermIterator<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "HpoParentIterator")
    }
}
