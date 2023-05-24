//! Utility structs and methods
use std::cmp::Ordering::{Equal, Less};

/// Iterator of all one-way pairwise combinations of the inner slice
///
/// # Examples
/// ```no_run
/// use hpo::utils::Combinations;
///
/// let mut c: Combinations<usize> = todo!(); // [1, 2, 3]
///
/// assert_eq!(c.next(), Some((&1, &2)));
/// assert_eq!(c.next(), Some((&1, &3)));
/// assert_eq!(c.next(), Some((&2, &3)));
/// assert!(c.next().is_none());
/// ```
pub struct Combinations<'a, T> {
    inner: &'a [Option<T>],
    idx1: usize,
    idx2: usize,
}

impl<'a, T> Combinations<'a, T> {
    /// Creates a new Combinations iterator
    pub(crate) fn new(inner: &'a [Option<T>]) -> Self {
        Self {
            inner,
            idx1: 0,
            idx2: 1,
        }
    }

    pub(crate) fn set_to_last(&mut self) {
        self.idx1 = self.inner.len() - 1;
        self.idx2 = 0;
    }
}

impl<'a, T> Iterator for Combinations<'a, T> {
    type Item = (&'a T, &'a T);
    fn next(&mut self) -> Option<Self::Item> {
        match (
            self.idx1 < self.inner.len(),
            self.idx2.cmp(&self.inner.len()),
        ) {
            (true, Less) => {
                self.idx2 += 1;
                if self.inner[self.idx1].is_none() {
                    return self.next();
                }
                if self.inner[self.idx2 - 1].is_none() {
                    return self.next();
                }
                Some((
                    self.inner[self.idx1]
                        .as_ref()
                        .expect("`None` values are skipped"),
                    self.inner[self.idx2 - 1]
                        .as_ref()
                        .expect("`None` values are skipped"),
                ))
            }
            (true, Equal) => {
                self.idx1 += 1;
                self.idx2 = self.idx1 + 1;
                self.next()
            }
            _ => None,
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn combinations() {
        let a = vec![Some(1), Some(2), Some(3), Some(4)];
        let mut c = Combinations::new(&a);
        assert_eq!(c.next(), Some((&1, &2)));
        assert_eq!(c.next(), Some((&1, &3)));
        assert_eq!(c.next(), Some((&1, &4)));
        assert_eq!(c.next(), Some((&2, &3)));
        assert_eq!(c.next(), Some((&2, &4)));
        assert_eq!(c.next(), Some((&3, &4)));
        assert_eq!(c.next(), None);
    }

    #[test]
    fn combinations_empty() {
        let a: Vec<Option<usize>> = vec![];
        let mut c = Combinations::new(&a);
        assert_eq!(c.next(), None);
    }

    #[test]
    fn combinations_single() {
        let a = vec![Some(1)];
        let mut c = Combinations::new(&a);
        assert_eq!(c.next(), None);
    }

    #[test]
    fn combinations_two() {
        let a = vec![Some(1), Some(2)];
        let mut c = Combinations::new(&a);
        assert_eq!(c.next(), Some((&1, &2)));
        assert_eq!(c.next(), None);
    }

    #[test]
    fn last_row() {
        let a = vec![Some(1), Some(2), Some(3), Some(4)];
        let mut c = Combinations::new(&a);
        c.set_to_last();
        assert_eq!(c.next(), Some((&4, &1)));
        assert_eq!(c.next(), Some((&4, &2)));
        assert_eq!(c.next(), Some((&4, &3)));
        assert_eq!(c.next(), Some((&4, &4)));
        assert_eq!(c.next(), None);
    }
}
