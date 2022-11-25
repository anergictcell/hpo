use std::collections::HashSet;
use std::ops::BitAnd;
use core::fmt::Debug;
use std::fmt::Display;


use crate::HpoError;
use crate::Ontology;
use crate::OntologyResult;

mod information_content;
pub use information_content::InformationContent;

mod internal;
pub use crate::term::internal::HpoTermInternal;

pub type HpoParents = HashSet<HpoTermId>;
pub type HpoChildren = HashSet<HpoTermId>;

#[derive(Copy, Clone, Eq, Hash, PartialEq)]
pub struct HpoTermId {
    inner: u32,
}

impl HpoTermId {
    fn new(s: &str) -> HpoTermId {
        HpoTermId {
            inner: s[3..].parse::<u32>().unwrap(),
        }
    }
}

impl From<&str> for HpoTermId {
    fn from(s: &str) -> Self {
        HpoTermId::new(s)
    }
}

impl From<usize> for HpoTermId {
    fn from(n: usize) -> Self {
        Self { inner: n as u32 }
    }
}

impl From<[char; 10]> for HpoTermId {
    fn from(s: [char; 10]) -> Self {
        let mut num = String::with_capacity(7);
        for c in &s[3..] {
            num.push(*c);
        }
        HpoTermId { inner: num.parse::<u32>().unwrap() }
    }
}

impl Debug for HpoTermId {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "HpoTermId({})", self)
    }
}

impl Display for HpoTermId {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.inner)
    }
}

impl PartialEq<str> for HpoTermId {
    fn eq(&self, other: &str) -> bool {
        self == &HpoTermId::new(other)
    }
}


#[derive(Debug)]
pub struct HpoTerm<'a> {
    id: &'a HpoTermId,
    name: &'a str,
    parents: &'a HpoParents,
    all_parents: &'a HpoParents,
    ontology: &'a Ontology,
}

impl<'a> HpoTerm<'a> {
    pub fn try_new(ontology: &'a Ontology, term: &HpoTermId) -> OntologyResult<HpoTerm<'a>> {
        let term = ontology.get(term).ok_or(HpoError::DoesNotExist)?;
        Ok(HpoTerm {
            id: term.id(),
            name: term.name(),
            parents: term.parents(),
            all_parents: term.all_parents(),
            ontology,
        })
    }

    pub fn new(ontology: &'a Ontology, term: &'a HpoTermInternal) -> HpoTerm<'a> {
        HpoTerm {
            id: term.id(),
            name: term.name(),
            parents: term.parents(),
            all_parents: term.all_parents(),
            ontology,
        }
    }

    pub fn id(&self) -> &HpoTermId {
        self.id
    }

    pub fn parents(&self) -> HpoParentIterator<'a> {
        HpoParentIterator::new(self.parents, self.ontology)
    }

    pub fn parent_ids(&self) -> &HpoParents {
        self.parents
    }

    pub fn all_parent_ids(&self) -> &HpoParents {
        self.all_parents
    }

    pub fn all_parents(&self) -> HpoParentIterator<'a> {
        HpoParentIterator::new(self.all_parents, self.ontology)
    }

    pub fn overlap(&self, other: &HpoTerm) -> HpoParents {
        self.all_parents.bitand(other.all_parents)
    }
}

pub struct HpoParentIterator<'a> {
    ontology: &'a Ontology,
    parents: std::collections::hash_set::Iter<'a, HpoTermId>,
}

impl<'a> HpoParentIterator<'a> {
    pub fn new(parents: &'a HpoParents, ontology: &'a Ontology) -> Self {
        HpoParentIterator {
            parents: parents.iter(),
            ontology,
        }
    }
}

impl<'a> std::iter::Iterator for HpoParentIterator<'a> {
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

impl Debug for HpoParentIterator<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "HpoParentIterator")
    }
}
