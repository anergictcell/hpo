use std::collections::HashSet;
use std::ops::{BitAnd, BitOr};

use crate::HpoTermId;

#[derive(Debug, Default, Clone)]
pub struct HpoGroup {
    ids: Vec<HpoTermId>,
}

impl HpoGroup {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            ids: Vec::with_capacity(capacity),
        }
    }

    pub fn is_empty(&self) -> bool {
        self.ids.is_empty()
    }

    pub fn len(&self) -> usize {
        self.ids.len()
    }

    pub fn insert(&mut self, id: HpoTermId) -> bool {
        match self.ids.binary_search(&id) {
            Ok(_) => false,
            Err(idx) => {
                self.ids.insert(idx, id);
                true
            }
        }
    }

    fn insert_unchecked(&mut self, id: HpoTermId) {
        self.ids.push(id)
    }

    pub fn contains(&self, id: &HpoTermId) -> bool {
        self.ids.binary_search(id).is_ok()
    }

    pub fn iter(&self) -> HpoGroupIterator {
        HpoGroupIterator::new(self.ids.iter())
    }

    // TODO: Is this really needed or makes sense?
    pub fn pop(&mut self) -> Option<HpoTermId> {
        self.ids.pop()
    }

    pub fn as_bytes(&self) -> Vec<u8> {
        let mut res = Vec::with_capacity(self.len() * 4);
        for id in &self.ids {
            res.append(&mut id.to_be_bytes().to_vec())
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

impl<'a> IntoIterator for &'a HpoGroup {
    type Item = &'a HpoTermId;

    type IntoIter = HpoGroupIterator<'a>;

    fn into_iter(self) -> HpoGroupIterator<'a> {
        HpoGroupIterator::new(self.ids.iter())
    }
}

pub struct HpoGroupIterator<'a> {
    inner: std::slice::Iter<'a, HpoTermId>,
}

impl<'a> HpoGroupIterator<'a> {
    fn new(inner: std::slice::Iter<'a, HpoTermId>) -> Self {
        Self { inner }
    }
}

impl<'a> Iterator for HpoGroupIterator<'a> {
    type Item = &'a HpoTermId;
    fn next(&mut self) -> Option<&'a HpoTermId> {
        self.inner.next()
    }
}

impl BitOr for &HpoGroup {
    type Output = HpoGroup;

    fn bitor(self, rhs: &HpoGroup) -> HpoGroup {
        let mut res = HpoGroup::with_capacity(self.len() + rhs.len());
        let (large, small) = if self.len() > rhs.len() {
            (self, rhs)
        } else {
            (rhs, self)
        };

        for id in &large.ids {
            res.insert_unchecked(*id);
        }
        for id in &small.ids {
            res.insert(*id);
        }
        res
    }
}

impl BitAnd for &HpoGroup {
    type Output = HpoGroup;

    fn bitand(self, rhs: &HpoGroup) -> HpoGroup {
        let mut res = HpoGroup::with_capacity(self.len());
        let (large, small) = if self.len() > rhs.len() {
            (self, rhs)
        } else {
            (rhs, self)
        };

        for id in &small.ids {
            if large.ids.contains(id) {
                res.insert_unchecked(*id);
            }
        }
        res
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
