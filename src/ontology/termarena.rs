use crate::term::internal::HpoTermInternal;
use crate::HpoTermId;
use std::{collections::{HashMap, hash_map::{Values, ValuesMut}}, hash::{BuildHasher, Hasher}};


struct HpoTermIdHasher {
    inner: u32
}

impl HpoTermIdHasher {
    fn new() -> Self {
        Self {inner: 0}
    }
}

impl Hasher for HpoTermIdHasher {
    fn finish(&self) -> u64 {
        self.inner.into()
    }

    fn write(&mut self, _bytes: &[u8]) {
        panic!("HpoTermIdHasher only supports u32");
    }

    fn write_u32(&mut self, inner: u32) {
        self.inner = inner
    }
}

struct HpoBuildHasher{}

impl HpoBuildHasher {
    fn new() -> Self {
        Self {}
    }
}

impl BuildHasher for HpoBuildHasher {
    type Hasher = HpoTermIdHasher;
    fn build_hasher(&self) -> Self::Hasher {
        HpoTermIdHasher::new()
    }
}

pub(crate) struct Arena {
    terms: HashMap<HpoTermId, HpoTermInternal, HpoBuildHasher>,
}

impl Arena {
    pub fn len(&self) -> usize {
        self.terms.len()
    }

    pub fn insert(&mut self, term: HpoTermInternal) {
        let id = term.id();
        self.terms.insert(*id, term);
    }

    pub fn get(&self, id: HpoTermId) -> Option<&HpoTermInternal> {
        self.terms.get(&id)
    }

    pub fn get_unchecked(&self, id: HpoTermId) -> &HpoTermInternal {
        self.terms.get(&id).unwrap()
    }

    pub fn get_mut(&mut self, id: HpoTermId) -> Option<&mut HpoTermInternal> {
        self.terms.get_mut(&id)
    }

    pub fn get_unchecked_mut(&mut self, id: HpoTermId) -> &mut HpoTermInternal {
        self.terms.get_mut(&id).unwrap()
    }

    pub fn values(&self) -> Values<'_, HpoTermId, HpoTermInternal> {
        self.terms.values()
    }

    pub fn values_mut(&mut self) -> ValuesMut<'_, HpoTermId, HpoTermInternal> {
        self.terms.values_mut()
    }

    pub fn keys(&self) -> Vec<HpoTermId> {
        self.terms.keys().copied().collect()
    }
}

impl Default for Arena {
    fn default() -> Self {
        Self {terms: HashMap::with_capacity_and_hasher(
            20_000,
            HpoBuildHasher::new()
        )}
    }
}
