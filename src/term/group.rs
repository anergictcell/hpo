use std::collections::HashSet;
use std::ops::{BitAnd, BitOr};

use crate::{HpoTerm, HpoTermId, Ontology};

/// A set of [`HpoTermId`] representing a group of HPO terms
///
/// Each term can occur only once in the group.
///
/// This group is used e.g. for having a set of parent or child HPO Terms
#[derive(Debug, Default, Clone)]
pub struct HpoGroup {
    ids: Vec<HpoTermId>,
}

impl HpoGroup {
    /// Constructs a new, empty [`HpoGroup`]
    pub fn new() -> Self {
        Self::default()
    }

    /// Constructs a new, empty [`HpoGroup`] with the given capacity
    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            ids: Vec::with_capacity(capacity),
        }
    }

    /// Returns `true` if the group contains no [`HpoTermId`]s
    pub fn is_empty(&self) -> bool {
        self.ids.is_empty()
    }

    /// Returns the number of [`HpoTermId`]s in the group
    pub fn len(&self) -> usize {
        self.ids.len()
    }

    /// Adds a new [`HpoTermId`] to the group
    ///
    /// Returns whether the `HpoTermId` was newly inserted. That is:
    ///
    /// - If the group did not previously contain this `HpoTermId`, true is returned.
    /// - If the group already contained this `HpoTermId`, false is returned.
    ///
    pub fn insert(&mut self, id: HpoTermId) -> bool {
        match self.ids.binary_search(&id) {
            Ok(_) => false,
            Err(idx) => {
                self.ids.insert(idx, id);
                true
            }
        }
    }

    /// Adds a new [`HpoTermId`] to the group
    ///
    /// # Note
    ///
    /// This method will not check if the `HpoTermId` already exists
    /// and will add it to the end of the vector. That means the internal
    /// sort order and uniqueness is not guaranteed.
    ///
    /// Using this method wrongly can have fatal effects on the correctness
    /// of the Ontology's functionality
    fn insert_unchecked(&mut self, id: HpoTermId) {
        self.ids.push(id);
    }

    /// Returns `true` if the group contains the [`HpoTermId`]
    pub fn contains(&self, id: &HpoTermId) -> bool {
        self.ids.binary_search(id).is_ok()
    }

    /// Returns an Iterator of the [`HpoTermId`]s inside the group
    pub fn iter(&self) -> HpoTermIds {
        HpoTermIds::new(self.ids.iter())
    }

    /// Returns the [`HpoTermId`] at the given index
    ///
    /// If the index is out of bounds, `None` is returned.
    fn get(&self, index: usize) -> Option<&HpoTermId> {
        self.ids.get(index)
    }

    /// Returns a byte representation of the group
    ///
    /// This method is only used for serializing the group
    /// to indicate term - parent relationships when creating
    /// the byte format of the Ontology
    pub fn as_bytes(&self) -> Vec<u8> {
        let mut res = Vec::with_capacity(self.len() * 4);
        for id in &self.ids {
            res.append(&mut id.to_be_bytes().to_vec());
        }
        res
    }
}

impl From<HashSet<HpoTermId>> for HpoGroup {
    fn from(s: HashSet<HpoTermId>) -> Self {
        let mut group = HpoGroup::with_capacity(s.len());
        for id in s {
            group.insert(id);
        }
        group
    }
}

/// Iterate [`HpoTerm`]s
///
/// This struct creates [`HpoTerm`]s from an owned [`HpoGroup`]
///
/// It is used in situations where we cannot reference the [`HpoGroup`]
/// because it is short-lived and must own it instead
pub struct GroupCombine<'a> {
    inner: HpoGroup,
    ontology: &'a Ontology,
    idx: usize,
}

impl<'a> GroupCombine<'a> {
    /// Constructs a new [`GroupCombine`] from an [`HpoGroup`] and a reference
    /// to the [`Ontology`]
    pub fn new(inner: HpoGroup, ontology: &'a Ontology) -> Self {
        Self {
            inner,
            idx: 0,
            ontology,
        }
    }
}

impl<'a> Iterator for GroupCombine<'a> {
    type Item = HpoTerm<'a>;
    fn next(&mut self) -> Option<Self::Item> {
        let index = self.idx;
        self.idx += 1;
        match self.inner.get(index) {
            Some(term_id) => self.ontology.hpo(*term_id),
            None => None,
        }
    }
}

impl<'a> IntoIterator for &'a HpoGroup {
    type Item = HpoTermId;

    type IntoIter = HpoTermIds<'a>;

    fn into_iter(self) -> HpoTermIds<'a> {
        HpoTermIds::new(self.ids.iter())
    }
}

/// An iterator over [`HpoTermId`]s
pub struct HpoTermIds<'a> {
    inner: std::slice::Iter<'a, HpoTermId>,
}

impl<'a> HpoTermIds<'a> {
    fn new(inner: std::slice::Iter<'a, HpoTermId>) -> Self {
        Self { inner }
    }
}

impl<'a> Iterator for HpoTermIds<'a> {
    type Item = HpoTermId;
    fn next(&mut self) -> Option<HpoTermId> {
        self.inner.next().copied()
    }
}

impl BitOr for &HpoGroup {
    type Output = HpoGroup;

    fn bitor(self, rhs: &HpoGroup) -> HpoGroup {
        let mut group = HpoGroup::with_capacity(self.len() + rhs.len());
        let (large, small) = if self.len() > rhs.len() {
            (self, rhs)
        } else {
            (rhs, self)
        };

        for id in &large.ids {
            group.insert_unchecked(*id);
        }
        for id in &small.ids {
            group.insert(*id);
        }
        group
    }
}

impl BitAnd for &HpoGroup {
    type Output = HpoGroup;

    fn bitand(self, rhs: &HpoGroup) -> HpoGroup {
        let mut group = HpoGroup::with_capacity(self.len());
        let (large, small) = if self.len() > rhs.len() {
            (self, rhs)
        } else {
            (rhs, self)
        };

        for id in &small.ids {
            if large.ids.contains(id) {
                group.insert_unchecked(*id);
            }
        }
        group
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hpogroup_iter() {
        let mut group = HpoGroup::new();
        group.insert(1u32.into());
        group.insert(2u32.into());
        group.insert(3u32.into());

        let mut ids = Vec::new();
        for id in &group {
            ids.push(id)
        }
        assert_eq!(ids.len(), 3);

        for id in &group {
            ids.push(id)
        }
        assert_eq!(ids.len(), 6);
    }

    #[test]
    fn test_bitor_set1() {
        let mut group1 = HpoGroup::new();
        group1.insert(1u32.into());
        group1.insert(2u32.into());
        group1.insert(3u32.into());

        let mut group2 = HpoGroup::new();
        group2.insert(2u32.into());
        group2.insert(4u32.into());

        let result = group1.bitor(&group2);
        let expected: Vec<HpoTermId> = vec![1u32.into(), 2u32.into(), 3u32.into(), 4u32.into()];
        assert_eq!(result.ids, expected);
    }

    #[test]
    fn test_bitor_set2() {
        let mut group1 = HpoGroup::new();
        group1.insert(1u32.into());
        group1.insert(2u32.into());
        group1.insert(3u32.into());

        let mut group2 = HpoGroup::new();
        group2.insert(1u32.into());
        group2.insert(2u32.into());
        group2.insert(4u32.into());
        group2.insert(5u32.into());

        let result = group1.bitor(&group2);
        let expected: Vec<HpoTermId> = vec![
            1u32.into(),
            2u32.into(),
            3u32.into(),
            4u32.into(),
            5u32.into(),
        ];
        assert_eq!(result.ids, expected);
    }

    #[test]
    fn test_bitand() {
        let mut group1 = HpoGroup::new();
        group1.insert(1u32.into());
        group1.insert(2u32.into());
        group1.insert(3u32.into());

        let mut group2 = HpoGroup::new();
        group2.insert(2u32.into());
        group2.insert(4u32.into());
        group2.insert(5u32.into());
        group2.insert(1u32.into());

        let result = group1.bitand(&group2);
        let expected: Vec<HpoTermId> = vec![1u32.into(), 2u32.into()];
        assert_eq!(result.ids, expected);
    }
}
