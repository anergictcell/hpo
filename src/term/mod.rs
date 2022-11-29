use std::collections::HashSet;
use core::fmt::Debug;


use crate::Ontology;

mod information_content;
mod internal;
mod parents;
mod hpoterm;
mod hpotermid;
mod group;

pub use information_content::InformationContent;
pub use internal::HpoTermInternal;
pub use hpoterm::HpoTerm;
pub use hpotermid::HpoTermId;
pub use group::{HpoGroup, HpoGroupIterator};

pub type HpoParents = HpoGroup;
pub type HpoChildren = HashSet<HpoTermId>;



pub struct HpoTermIterator<'a> {
    ontology: &'a Ontology,
    parents: HpoGroupIterator<'a>,
}

impl<'a> HpoTermIterator<'a> {
    pub fn new(parents: &'a HpoParents, ontology: &'a Ontology) -> Self {
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
                let term = self.ontology.get(term).unwrap();
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
