//! This module contains the [`HpoGroup`] struct and the corresponding [`HpoTermId`] iterators [`Iter`] and [`Combined`]
use crate::annotations::AnnotationId;
use std::collections::HashSet;
use std::ops::{Add, BitAnd, BitOr};

use smallvec::SmallVec;

use crate::term;
use crate::{HpoTerm, HpoTermId, Ontology};

/// A set of [`HpoTermId`] representing a group of HPO terms
///
/// Each term can occur only once in the group.
///
/// This group is used e.g. for having a set of parent or child HPO Terms
#[derive(Debug, Default, Clone)]
pub struct HpoGroup {
    ids: SmallVec<[HpoTermId; 30]>,
}

impl HpoGroup {
    /// Constructs a new, empty [`HpoGroup`]
    pub fn new() -> Self {
        Self::default()
    }

    /// Constructs a new, empty [`HpoGroup`] with the given capacity
    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            ids: SmallVec::with_capacity(capacity),
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

    /// Removes all [`HpoTermId`] from the group and empties it
    pub fn clear(&mut self) {
        self.ids.clear();
    }

    /// Returns an Iterator of the [`HpoTermId`]s inside the group
    pub fn iter(&self) -> Iter {
        self.into_iter()
    }

    /// Returns the [`HpoTermId`] at the given index
    ///
    /// If the index is out of bounds, `None` is returned.
    pub fn get(&self, index: usize) -> Option<&HpoTermId> {
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

    /// Returns an iterator of [`HpoTerm`]
    pub fn terms<'a>(&'a self, ontology: &'a Ontology) -> term::Iter<'a> {
        term::Iter::new(self.iter(), ontology)
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

impl From<Vec<HpoTermId>> for HpoGroup {
    fn from(s: Vec<HpoTermId>) -> Self {
        let mut group = HpoGroup::with_capacity(s.len());
        for id in s {
            group.insert(id);
        }
        group
    }
}

impl From<Vec<u32>> for HpoGroup {
    fn from(s: Vec<u32>) -> Self {
        let mut group = HpoGroup::with_capacity(s.len());
        for id in s {
            group.insert(id.into());
        }
        group
    }
}

impl FromIterator<HpoTermId> for HpoGroup {
    fn from_iter<T: IntoIterator<Item = HpoTermId>>(iter: T) -> Self {
        let mut group = HpoGroup::new();
        for id in iter {
            group.insert(id);
        }
        group
    }
}

impl<'a> FromIterator<HpoTerm<'a>> for HpoGroup {
    fn from_iter<T: IntoIterator<Item = HpoTerm<'a>>>(iter: T) -> Self {
        let mut group = HpoGroup::new();
        for term in iter {
            group.insert(term.id());
        }
        group
    }
}

impl<'a> IntoIterator for &'a HpoGroup {
    type Item = HpoTermId;

    type IntoIter = Iter<'a>;
    fn into_iter(self) -> Self::IntoIter {
        Iter {
            iter: self.ids.iter(),
        }
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

impl BitOr for HpoGroup {
    type Output = HpoGroup;

    fn bitor(self, rhs: HpoGroup) -> HpoGroup {
        (&self).bitor(&rhs)
    }
}

impl BitOr<&HpoGroup> for HpoGroup {
    type Output = HpoGroup;

    fn bitor(self, rhs: &HpoGroup) -> HpoGroup {
        (&self).bitor(rhs)
    }
}

impl BitOr<HpoTermId> for &HpoGroup {
    type Output = HpoGroup;

    fn bitor(self, rhs: HpoTermId) -> HpoGroup {
        let mut group = self.clone();
        group.insert(rhs);
        group
    }
}

impl Add<HpoTermId> for &HpoGroup {
    type Output = HpoGroup;
    fn add(self, rhs: HpoTermId) -> Self::Output {
        let mut group = self.clone();
        group.insert(rhs);
        group
    }
}

impl Add<HpoTermId> for HpoGroup {
    type Output = HpoGroup;
    fn add(self, rhs: HpoTermId) -> Self::Output {
        (&self) + rhs
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

impl BitAnd for HpoGroup {
    type Output = HpoGroup;

    fn bitand(self, rhs: HpoGroup) -> HpoGroup {
        (&self).bitand(&rhs)
    }
}

impl BitAnd<&HpoGroup> for HpoGroup {
    type Output = HpoGroup;

    fn bitand(self, rhs: &HpoGroup) -> HpoGroup {
        (&self).bitand(rhs)
    }
}

/// [`HpoTermId`] iterator
pub struct Iter<'a> {
    iter: std::slice::Iter<'a, HpoTermId>,
}

impl Iterator for Iter<'_> {
    type Item = HpoTermId;
    fn next(&mut self) -> Option<Self::Item> {
        self.iter.next().copied()
    }
}

/// [`HpoTerm`] iterator for an owned [`HpoGroup`]
///
/// This iterator is needed for some cases where an `HpoGroup` is created
/// by a method and does not live long enough to be used with [`Iter`].
pub struct Combined<'a> {
    group: HpoGroup,
    ontology: &'a Ontology,
}

impl<'a> Combined<'a> {
    /// Constructs a new [`Combined`] from an [`HpoGroup`] and a reference
    /// to the [`Ontology`]
    pub fn new(group: HpoGroup, ontology: &'a Ontology) -> Self {
        Self { group, ontology }
    }

    /// Returns an Iterator of [`HpoTermId`]
    pub fn iter(&self) -> term::Iter<'_> {
        self.into_iter()
    }

    /// Returns the number of [`HpoTermId`] in the Iterator
    pub fn len(&self) -> usize {
        self.group.len()
    }

    /// Returns true if the Iterator does not contain any [`HpoTermId`]
    pub fn is_empty(&self) -> bool {
        self.group.is_empty()
    }
}

impl<'a> IntoIterator for &'a Combined<'a> {
    /// iterates [`HpoTerm`]s
    type Item = HpoTerm<'a>;

    /// [`HpoTerm`] Iterator
    type IntoIter = term::Iter<'a>;

    fn into_iter(self) -> Self::IntoIter {
        term::Iter::new(self.group.iter(), self.ontology)
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
        assert_eq!(result.ids.into_vec(), expected);
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
        assert_eq!(result.ids.into_vec(), expected);
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

        let result = &group1 & &group2;
        let expected: Vec<HpoTermId> = vec![1u32.into(), 2u32.into()];
        assert_eq!(result.ids.into_vec(), expected);
    }
}
