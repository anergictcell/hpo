use crate::HpoTerm;
use crate::term::HpoTermIterator;
use crate::Ontology;

pub struct SetSimilarity {}

pub struct HpoSet<'a> {
    inner: Vec<HpoTerm<'a>>
}

impl<'a> HpoSet<'a> {
    pub fn child_nodes(&mut self) -> HpoSet {
        self.inner.iter().filter(|term1| {
            let mut has_children = false;
            for term2 in &self.inner {
                if term1.parent_of(term2) {
                    has_children = true
                }
            }
            !has_children
        }).copied().collect::<Vec<HpoTerm<'a>>>().into()
    }

    pub fn len(&self) -> usize {
        self.inner.len()
    }

    pub fn remove_modifier(&mut self) {
        unimplemented!()
    }

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


impl<'a> From<HpoTermIterator<'a>> for HpoSet<'a> {
    fn from(iter: HpoTermIterator<'a>) -> HpoSet<'a> {
        Self {inner: iter.collect()}
    }
}


impl<'a> From<Vec<HpoTerm<'a>>> for HpoSet<'a> {
    fn from(terms: Vec<HpoTerm<'a>>) -> Self {
        Self {inner: terms}
    }
}