//! A custom matrix for quick row and column-based data access
//!
//! `Matrix` is only used within the [`SimilarityCombiner`](`crate::similarity::SimilarityCombiner`) trait
//! and should not be used for other purposes, as it does not contain many safety guarantees.
//!
//! Imagine the following matrix / Dataframe
//!
//! | Index |   0 |   1 |   2 |   3 |
//! |:----- | ---:| ---:| ---:| ---:|
//! | **0** |  11 |  12 |  13 |  14 |
//! | **1** |  21 |  22 |  23 |  24 |
//! | **2** |  31 |  32 |  33 |  34 |
//! 
//! ```no_run
//! use hpo::matrix::Matrix;
//! let data = vec![11, 12, 13, 14, 21, 22, 23, 24, 31, 32, 33, 34];
//! let m = Matrix::new(3, 4, &data);
//!
//! for row in m.rows() {
//!     let v: Vec<String> = row.map(|v| format!("{}", v)).collect();
//!     println!("[{}]", v.join(", "));
//! }
//!
//! // >> [11, 12, 13, 14]
//! // >> [21, 22, 23, 24]
//! // >> [31, 32, 33, 34]
//!
//!
//! for col in m.cols() {
//!     let v: Vec<String> = col.map(|v| format!("{}", v)).collect();
//!     println!("[{}]", v.join(", "));
//! }
//!
//! // >> [11, 21, 31]
//! // >> [12, 22, 32]
//! // >> [13, 23, 33]
//! // >> [14, 24, 34]
//! ```
//!
//! There are no logic checks to ensure that the rows and column
//! match the data length, so callers must ensure this
//!
use std::fmt::Debug;

/// A custom matrix for quick row and column-based data access
///
/// This struct is used within the [`SimilarityCombiner`](`crate::similarity::SimilarityCombiner`) trait
pub struct Matrix<'a, T> {
    rows: usize,
    cols: usize,
    data: &'a [T],
}

impl<'a, T> Matrix<'a, T> {
    /// Crates a new Matrix from `data` with `rows` number of rows, `cols` number of columns
    ///
    /// # Note:
    /// It does not check if the number of rows and cols correspond to the length of data.
    ///
    /// # Examples
    /// ```rust
    /// use hpo::matrix::Matrix;
    /// let data = vec![11, 12, 13, 14, 21, 22, 23, 24, 31, 32, 33, 34];
    /// let m = Matrix::new(3, 4, &data);
    /// assert_eq!(m.len(), 12);
    /// ```
    pub fn new(rows: usize, cols: usize, data: &'a [T]) -> Self {
        Self { rows, cols, data }
    }

    /// Returns the total length of the data
    ///
    /// # Examples
    /// ```rust
    /// use hpo::matrix::Matrix;
    /// let data = vec![11, 12, 13, 14, 21, 22, 23, 24, 31, 32, 33, 34];
    /// let m = Matrix::new(3, 4, &data);
    /// assert_eq!(m.len(), 12);
    /// ```
    pub fn len(&self) -> usize {
        self.data.len()
    }

    /// Returns true if the `Matrix` does not contain any data
    ///
    /// # Examples
    /// ```rust
    /// use hpo::matrix::Matrix;
    /// let data: Vec<usize> = vec![];
    /// let m = Matrix::new(0, 0, &data);
    /// assert!(m.is_empty());
    /// ```
    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    /// Returns a Tuple with number of rows and number of columns
    ///
    /// # Examples
    /// ```rust
    /// use hpo::matrix::Matrix;
    /// let data = vec![11, 12, 13, 14, 21, 22, 23, 24, 31, 32, 33, 34];
    /// let m = Matrix::new(3, 4, &data);
    /// assert_eq!(m.dim(), (3, 4));
    /// ```
    pub fn dim(&self) -> (usize, usize) {
        (self.rows, self.cols)
    }

    /// Iterates the rows of the matrix
    ///
    /// # Examples
    /// ```rust
    /// use hpo::matrix::Matrix;
    /// let data = vec![11, 12, 13, 14, 21, 22, 23, 24, 31, 32, 33, 34];
    /// let m = Matrix::new(3, 4, &data);
    ///
    /// let mut rows = m.rows();
    ///
    /// let mut row = rows.next().unwrap();
    /// assert_eq!(row.next().unwrap(), &11);
    /// assert_eq!(row.next().unwrap(), &12);
    /// assert_eq!(row.next().unwrap(), &13);
    /// assert_eq!(row.next().unwrap(), &14);
    /// assert!(row.next().is_none());
    ///
    /// assert_eq!(rows.next().unwrap().count(), 4);
    /// assert_eq!(rows.next().unwrap().count(), 4);
    /// assert!(rows.next().is_none());
    /// ```
    pub fn rows(&self) -> RowIterator<T> {
        RowIterator::new(self.data, self.row_indicies())
    }

    /// Iterates the columns of the matrix
    ///
    /// # Examples
    /// ```rust
    /// use hpo::matrix::Matrix;
    /// let data = vec![11, 12, 13, 14, 21, 22, 23, 24, 31, 32, 33, 34];
    /// let m = Matrix::new(3, 4, &data);
    ///
    /// let mut cols = m.cols();
    ///
    /// let mut col = cols.next().unwrap();
    /// assert_eq!(col.next().unwrap(), &11);
    /// assert_eq!(col.next().unwrap(), &21);
    /// assert_eq!(col.next().unwrap(), &31);
    /// assert!(col.next().is_none());
    ///
    /// assert_eq!(cols.next().unwrap().count(), 3);
    /// assert_eq!(cols.next().unwrap().count(), 3);
    /// assert_eq!(cols.next().unwrap().count(), 3);
    /// assert!(cols.next().is_none());
    /// ```
    pub fn cols(&self) -> ColumnIterator<T> {
        ColumnIterator::new(self.data, self.col_indicies())
    }

    fn row_indicies(&self) -> RowIndexIterator {
        RowIndexIterator::new(self.rows, self.cols)
    }

    fn col_indicies(&self) -> ColumnIndexIterator {
        ColumnIndexIterator::new(self.cols)
    }
}

impl<T: std::fmt::Display> Debug for Matrix<'_, T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for row in self.rows() {
            let v: Vec<String> = row.map(|v| format!("{v}")).collect();
            writeln!(f, "[{}]", v.join(", ")).unwrap();
        }
        Ok(())
    }
}

/// An iterator of the values of a single row of a `Matrix`
///
/// This struct is yielded by `RowIterator`
pub struct Row<'a, T> {
    iter: std::slice::Iter<'a, T>,
}

impl<'a, T> Row<'a, T> {
    fn new(iter: std::slice::Iter<'a, T>) -> Self {
        Row { iter }
    }
}

impl<'a, T> Iterator for Row<'a, T> {
    type Item = &'a T;
    fn next(&mut self) -> Option<Self::Item> {
        self.iter.next()
    }
}

/// Iterates the rows of a `Matrix`, returning an Iterator over individual row values
///
/// This struct is yielded by `Matrix::rows`
pub struct RowIterator<'a, T> {
    iter: RowIndexIterator,
    data: &'a [T],
}

impl<'a, T> RowIterator<'a, T> {
    fn new(data: &'a [T], iter: RowIndexIterator) -> Self {
        Self { iter, data }
    }
}

impl<'a, T> Iterator for RowIterator<'a, T> {
    type Item = Row<'a, T>;
    fn next(&mut self) -> Option<Self::Item> {
        match self.iter.next() {
            Some(range) => Some(Row::new(self.data[range].iter())),
            None => None,
        }
    }
}

/// Yields a tuple with start and end position for each row
struct RowIndexIterator {
    rows: usize,
    cols: usize,
    idx: usize,
}

impl RowIndexIterator {
    fn new(rows: usize, cols: usize) -> Self {
        Self { rows, cols, idx: 0 }
    }
}

impl Iterator for RowIndexIterator {
    type Item = std::ops::RangeInclusive<usize>;
    fn next(&mut self) -> Option<Self::Item> {
        if self.idx >= self.rows * self.cols {
            return None;
        }
        let res = std::ops::RangeInclusive::new(self.idx, self.idx + self.cols - 1);
        self.idx += self.cols;
        Some(res)
    }
}

/// An iterator of the values of a single column of a `Matrix`
///
/// This struct is yielded by `ColumnIterator`
pub struct Column<'a, T> {
    iter: std::iter::StepBy<std::slice::Iter<'a, T>>,
}

impl<'a, T> Column<'a, T> {
    fn new(iter: std::iter::StepBy<std::slice::Iter<'a, T>>) -> Self {
        Self { iter }
    }
}

impl<'a, T> Iterator for Column<'a, T> {
    type Item = &'a T;
    fn next(&mut self) -> Option<Self::Item> {
        self.iter.next()
    }
}

/// Iterates the columns of a `Matrix`, returning an Iterator over individual column values
///
/// This struct is yielded by `Matrix::cols`
pub struct ColumnIterator<'a, T> {
    data: &'a [T],
    iter: ColumnIndexIterator,
}

impl<'a, T> ColumnIterator<'a, T> {
    fn new(data: &'a [T], iter: ColumnIndexIterator) -> Self {
        Self { data, iter }
    }
}

impl<'a, T> Iterator for ColumnIterator<'a, T> {
    type Item = Column<'a, T>;
    fn next(&mut self) -> Option<Self::Item> {
        match self.iter.next() {
            None => None,
            Some(range) => {
                let mut iter = self.data.iter();
                // iter.advance_by(range.idx); // Not yet stable
                for _ in 0..range.idx {
                    iter.next();
                }
                Some(Column::new(iter.step_by(range.step)))
            }
        }
    }
}

#[derive(Debug, PartialEq, Eq)]
struct ColumnRange {
    idx: usize,
    step: usize,
}

impl ColumnRange {
    pub fn new(idx: usize, step: usize) -> Self {
        Self { idx, step }
    }
}

struct ColumnIndexIterator {
    cols: usize,
    idx: usize,
}

impl ColumnIndexIterator {
    fn new(cols: usize) -> Self {
        Self { cols, idx: 0 }
    }
}

impl Iterator for ColumnIndexIterator {
    type Item = ColumnRange;
    fn next(&mut self) -> Option<Self::Item> {
        if self.idx >= self.cols {
            return None;
        }
        let res = ColumnRange::new(self.idx, self.cols);
        self.idx += 1;
        Some(res)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_row_index_iterator() {
        let mut iter = RowIndexIterator::new(2, 4);
        assert_eq!(iter.next(), Some(0..=3));
        assert_eq!(iter.next(), Some(4..=7));
        assert_eq!(iter.next(), None);

        let mut iter = RowIndexIterator::new(3, 2);
        assert_eq!(iter.next(), Some(0..=1));
        assert_eq!(iter.next(), Some(2..=3));
        assert_eq!(iter.next(), Some(4..=5));
        assert_eq!(iter.next(), None);

        let mut iter = RowIndexIterator::new(1, 10);
        assert_eq!(iter.next(), Some(0..=9));
        assert_eq!(iter.next(), None);

        let mut iter = RowIndexIterator::new(0, 10);
        assert_eq!(iter.next(), None);

        let mut iter = RowIndexIterator::new(2, 4);
        let mut row = iter.next().unwrap();
        // row == 0..=3
        assert_eq!(row.next(), Some(0));
        assert_eq!(row.next(), Some(1));
        assert_eq!(row.next(), Some(2));
        assert_eq!(row.next(), Some(3));
        assert_eq!(row.next(), None);
    }

    #[test]
    fn test_row_iterator() {
        let data = vec![1, 2, 3, 4, 5, 6];
        let mut rowiter = RowIterator::new(&data, RowIndexIterator::new(2, 3));

        let mut row = rowiter.next().unwrap();
        assert_eq!(row.next(), Some(&1));
        assert_eq!(row.next(), Some(&2));
        assert_eq!(row.next(), Some(&3));
        assert_eq!(row.next(), None);

        let mut row = rowiter.next().unwrap();
        assert_eq!(row.next(), Some(&4));
        assert_eq!(row.next(), Some(&5));
        assert_eq!(row.next(), Some(&6));
        assert_eq!(row.next(), None);

        assert!(rowiter.next().is_none());
    }

    #[test]
    fn test_empty_row_iterator() {
        let data = vec![1, 2, 3, 4, 5, 6];
        let mut rowiter = RowIterator::new(&data, RowIndexIterator::new(0, 3));
        assert!(rowiter.next().is_none());
    }

    #[test]
    fn test_row_iterator_sums() {
        let data = vec![1, 2, 3, 4, 5, 6];
        let mut rowiter = RowIterator::new(&data, RowIndexIterator::new(2, 3));

        let row = rowiter.next().unwrap();
        assert_eq!(row.sum::<i32>(), 6);

        let row = rowiter.next().unwrap();
        assert_eq!(row.sum::<i32>(), 15);
    }

    #[test]
    fn test_column_index_iterator() {
        let mut iter = ColumnIndexIterator::new(4);
        assert_eq!(iter.next(), Some(ColumnRange::new(0, 4)));
        assert_eq!(iter.next(), Some(ColumnRange::new(1, 4)));
        assert_eq!(iter.next(), Some(ColumnRange::new(2, 4)));
        assert_eq!(iter.next(), Some(ColumnRange::new(3, 4)));
        assert_eq!(iter.next(), None);

        let mut iter = ColumnIndexIterator::new(2);
        assert_eq!(iter.next(), Some(ColumnRange::new(0, 2)));
        assert_eq!(iter.next(), Some(ColumnRange::new(1, 2)));
        assert_eq!(iter.next(), None);

        let mut iter = ColumnIndexIterator::new(1);
        assert_eq!(iter.next(), Some(ColumnRange::new(0, 1)));
        assert_eq!(iter.next(), None);

        let mut iter = ColumnIndexIterator::new(0);
        assert_eq!(iter.next(), None);
    }

    #[test]
    fn test_column_iterator() {
        let data = vec![1, 2, 3, 4, 5, 6];
        let mut coliter = ColumnIterator::new(&data, ColumnIndexIterator::new(3));

        let mut col = coliter.next().unwrap();
        assert_eq!(col.next(), Some(&1));
        assert_eq!(col.next(), Some(&4));
        assert_eq!(col.next(), None);

        let mut col = coliter.next().unwrap();
        assert_eq!(col.next(), Some(&2));
        assert_eq!(col.next(), Some(&5));
        assert_eq!(col.next(), None);

        let mut col = coliter.next().unwrap();
        assert_eq!(col.next(), Some(&3));
        assert_eq!(col.next(), Some(&6));
        assert_eq!(col.next(), None);

        assert!(coliter.next().is_none());
    }

    #[test]
    fn test_col_iterator_sums() {
        let data = vec![1, 2, 3, 4, 5, 6];
        let mut coliter = ColumnIterator::new(&data, ColumnIndexIterator::new(3));

        let col = coliter.next().unwrap();
        assert_eq!(col.sum::<i32>(), 5);

        let col = coliter.next().unwrap();
        assert_eq!(col.sum::<i32>(), 7);

        let col = coliter.next().unwrap();
        assert_eq!(col.sum::<i32>(), 9);
    }
}
