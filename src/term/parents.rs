use crate::HpoTermId;

const TERM_ARRAY_SIZE: usize = 10_000_000;

struct HpoParents {
    inner: [bool; TERM_ARRAY_SIZE]
}

impl Default for HpoParents {
    fn default() -> Self {
        Self {inner: [false; TERM_ARRAY_SIZE]}
    }
}

impl IntoIterator for HpoParents {
    type Item = HpoTermId;
    type IntoIter = HpoParentIdIterator;
    fn into_iter(self) -> HpoParentIdIterator {
        HpoParentIdIterator::new(&self.inner)
    }
}


struct HpoParentIdIterator {
    idx: usize,
    parents: [bool; TERM_ARRAY_SIZE]
}


impl HpoParentIdIterator {
    pub fn new(parents: &[bool; TERM_ARRAY_SIZE]) -> Self {
        Self { idx: 0, parents: *parents }
    }
}

impl Iterator for HpoParentIdIterator {
    type Item = HpoTermId;
    fn next(&mut self) -> Option<Self::Item> {
        loop {
            self.idx += 1;
            if self.idx >= TERM_ARRAY_SIZE {
                return None;
            }
            if self.parents[self.idx] {
                return Some(HpoTermId::from(self.idx));
            }
        }
       
    }
}