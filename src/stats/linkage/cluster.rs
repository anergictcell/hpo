//! A Cluster of `HpoSet` for linkage analysis.
//!
//! A cluster contains at least 2 `HpoSet`, but can also group other cluster.

/// A combination of (at least) 2 `HpoSet` for Dendogram clustering
///
/// This `struct` is created by the [`Linkage`](`crate::stats::Linkage`) struct and yielded
/// by the [`Linkage::into_cluster`](`crate::stats::Linkage::into_cluster`) and 
/// [`Linkage::cluster`](`crate::stats::Linkage::cluster`) iterators
#[derive(Debug, Clone, Copy)]
pub struct Cluster {
    idx1: usize,
    idx2: usize,
    distance: f32,
    size: usize,
}

impl Cluster {
    /// Creates a new `Cluster`
    pub(super) fn new(idx1: usize, idx2: usize, distance: f32, size: usize) -> Self {
        Self {
            idx1,
            idx2,
            distance,
            size,
        }
    }

    /// Returns the index of the left hand side node
    pub fn lhs(&self) -> usize {
        self.idx1
    }

    /// Returns the index of the right hand side node
    pub fn rhs(&self) -> usize {
        self.idx2
    }

    /// Returns the distance of the two nodes
    pub fn distance(&self) -> f32 {
        self.distance
    }

    /// Returns the total number of child nodes
    #[allow(clippy::len_without_is_empty)]
    pub fn len(&self) -> usize {
        self.size
    }
}

#[derive(Debug, Default)]
pub(super) struct ClusterVec(Vec<Cluster>);

impl ClusterVec {
    pub(super) fn with_capacity(capacity: usize) -> Self {
        Self(Vec::with_capacity(capacity))
    }

    pub(super) fn push(&mut self, value: Cluster) {
        self.0.push(value);
    }

    pub(super) fn get(&self, index: usize) -> Option<&Cluster> {
        self.0.get(index)
    }

    pub(super) fn iter(&self) -> Iter<'_> {
        self.into_iter()
    }
}

impl IntoIterator for ClusterVec {
    type Item = Cluster;
    type IntoIter = IntoIter;
    fn into_iter(self) -> Self::IntoIter {
        IntoIter::new(self.0.into_iter())
    }
}

impl<'a> IntoIterator for &'a ClusterVec {
    type Item = &'a Cluster;
    type IntoIter = Iter<'a>;
    fn into_iter(self) -> Self::IntoIter {
        Iter::new(self.0.iter())
    }
}

/// Iterates [`Cluster`] references
///
/// This `struct` is yielded by the [`Linkage::cluster`](`crate::stats::Linkage::cluster`) method
pub struct Iter<'a> {
    iter: std::slice::Iter<'a, Cluster>,
}

impl<'a> Iter<'a> {
    fn new(iter: std::slice::Iter<'a, Cluster>) -> Self {
        Self { iter }
    }
}

impl<'a> Iterator for Iter<'a> {
    type Item = &'a Cluster;

    fn next(&mut self) -> Option<Self::Item> {
        self.iter.next()
    }
}

impl DoubleEndedIterator for Iter<'_> {
    fn next_back(&mut self) -> Option<Self::Item> {
        self.iter.next_back()
    }
}

impl ExactSizeIterator for Iter<'_> {
    fn len(&self) -> usize {
        self.iter.len()
    }
}


/// Iterates [`Cluster`]
///
/// This `struct` is yielded by the [`Linkage::into_cluster`](`crate::stats::Linkage::into_cluster`) method
#[derive(Debug)]
pub struct IntoIter {
    iter: std::vec::IntoIter<Cluster>,
}

impl IntoIter {
    fn new(iter: std::vec::IntoIter<Cluster>) -> Self {
        Self { iter }
    }
}

impl Iterator for IntoIter {
    type Item = Cluster;
    fn next(&mut self) -> Option<Self::Item> {
        self.iter.next()
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        self.iter.size_hint()
    }
}

impl DoubleEndedIterator for IntoIter {
    fn next_back(&mut self) -> Option<Self::Item> {
        self.iter.next_back()
    }
}

impl ExactSizeIterator for IntoIter {
    fn len(&self) -> usize {
        self.iter.len()
    }
}
