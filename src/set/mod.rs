//! An HpoSet can represent e.g. the clinical information of a patient or the symptoms of a disease
use crate::term::HpoTerms;
use crate::HpoTerm;
use crate::Ontology;

/// The HpoSet contains a set of unique HPO terms
///
/// As in a set, each term can only appear once
/// though that is not yet guaranteed in the implementation (TODO)
pub struct HpoSet<'a> {
    inner: Vec<HpoTerm<'a>>,
}

impl<'a> HpoSet<'a> {
    /// Returns a new HpoSet that contains only the child-most terms
    ///
    /// This means that it only contains terms that don't have a child
    /// term present in the set.
    pub fn child_nodes(&mut self) -> HpoSet {
        self.inner
            .iter()
            .filter(|term1| {
                let mut has_children = false;
                for term2 in &self.inner {
                    if term1.parent_of(term2) {
                        has_children = true
                    }
                }
                !has_children
            })
            .copied()
            .collect::<Vec<HpoTerm<'a>>>()
            .into()
    }

    /// Returns the number of terms in the set
    pub fn len(&self) -> usize {
        self.inner.len()
    }

    /// Returns a new set without modifier terms
    ///
    /// This is not yet implemented
    pub fn remove_modifier(&mut self) {
        unimplemented!()
    }

    /// Returns a new set where all obsolete terms are replaced
    ///
    /// This is not yet implemented
    pub fn replace_obsolete(&mut self, _ontology: &Ontology) {
        unimplemented!()
    }
}

impl<'a> IntoIterator for &'a HpoSet<'a> {
    type Item = &'a HpoTerm<'a>;
    type IntoIter = std::slice::Iter<'a, HpoTerm<'a>>;
    fn into_iter(self) -> Self::IntoIter {
        self.inner.iter()
    }
}

impl<'a> From<HpoTerms<'a>> for HpoSet<'a> {
    fn from(iter: HpoTerms<'a>) -> HpoSet<'a> {
        Self {
            inner: iter.collect(),
        }
    }
}

impl<'a> From<Vec<HpoTerm<'a>>> for HpoSet<'a> {
    fn from(terms: Vec<HpoTerm<'a>>) -> Self {
        Self { inner: terms }
    }
}
