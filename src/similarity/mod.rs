//! Methods to calculate the Similarity between two terms or sets of terms

use std::cell::RefCell;
use std::collections::HashMap;

use crate::matrix::Matrix;
use crate::set::HpoSet;
use crate::term::InformationContentKind;
use crate::{HpoTerm, HpoTermId};

mod defaults;
pub use defaults::{Distance, GraphIc, InformationCoefficient, Jc, Lin, Relevance, Resnik};

/// Trait for similarity score calculation between 2 [`HpoTerm`]s
///
/// `hpo` already comes pre-loaded with several common and well established
/// similarity algorithms that implement the `Similarity` trait.
pub trait Similarity {
    /// calculates the actual similarity between term a and term b
    fn calculate(&self, a: &HpoTerm, b: &HpoTerm) -> f32;
}

/// This trait is needed for custom implementations
///
/// For similarity calculation between sets of `HpoTerm`s
/// the similarity scores must be combined
pub trait SimilarityCombiner {
    /// This method implements the actual logic to calculate a single
    /// similarity score from a Matrix of term - term similarity scores.
    fn combine(&self, m: &Matrix<f32>) -> f32;

    /// this method is called by [`GroupSimilarity`] to combine individual term - term
    /// similarity scores into a single score for the group - group similarity
    fn calculate(&self, m: &Matrix<f32>) -> f32 {
        if m.is_empty() {
            return 0.0;
        }
        self.combine(m)
    }

    /// Returns the maximum values of each row
    fn row_maxes(&self, m: &Matrix<f32>) -> Vec<f32> {
        m.rows()
            .map(|row| {
                // I have no idea why, but I could not get a.max(b) to work
                // with the borrow checker
                row.reduce(|a, b| if a > b { a } else { b }).unwrap()
            })
            .copied()
            .collect()
    }

    /// Returns the maximum values of each column
    fn col_maxes(&self, m: &Matrix<f32>) -> Vec<f32> {
        m.cols()
            .map(|col| {
                // I have no idea why, but I could not get a.max(b) to work
                // with the borrow checker
                col.reduce(|a, b| if a > b { a } else { b }).unwrap()
            })
            .copied()
            .collect()
    }

    /// Returns the dimenension of the `Matrix`, (rows, columns)
    fn dim_f32(&self, m: &Matrix<f32>) -> (f32, f32) {
        let (rows, cols) = m.dim();
        (usize_to_f32(rows), usize_to_f32(cols))
    }
}

/// Caches the Similarity score for each [`HpoTerm`] pair
///
/// Use this struct to wrap your Similarity method if you are
/// running many batch comparisons where it's highly likely that
/// several comparisons will be repeatedly run.
/// This is very useful when you're e.g. comparing the set of a patient
/// with every disease or gene.
///
/// # Note
///
/// This struct cannot be used in multithreaded processing
pub struct CachedSimilarity<T> {
    similarity: T,
    cache: RefCell<HashMap<(HpoTermId, HpoTermId), f32>>,
}

impl<T: Similarity> CachedSimilarity<T> {
    /// Constructs a new [`CachedSimilarity`] struct
    pub fn new(similarity: T) -> Self {
        Self {
            similarity,
            cache: RefCell::new(HashMap::default()),
        }
    }
}

impl<T: Similarity> Similarity for CachedSimilarity<T> {
    fn calculate(&self, a: &HpoTerm, b: &HpoTerm) -> f32 {
        *self
            .cache
            .borrow_mut()
            .entry((a.id(), b.id()))
            .or_insert_with(|| self.similarity.calculate(a, b))
    }
}

/// Default implementations for combining similarity scores
/// for comparison of 2 sets of terms
pub enum StandardCombiner {
    /// funSimAvg algorithm from [Schlicker A, et. al., BMC Bioinf (2006)](https://pubmed.ncbi.nlm.nih.gov/16776819/)
    FunSimAvg,
    /// funSimMax algorithm from [Schlicker A, et. al., BMC Bioinf (2006)](https://pubmed.ncbi.nlm.nih.gov/16776819/)
    FunSimMax,
    /// BMA algorithm from [Wang JZ, et. al., Bioinformatics (2007)](https://pubmed.ncbi.nlm.nih.gov/17344234/)
    Bwa,
}

impl Default for StandardCombiner {
    fn default() -> Self {
        Self::FunSimAvg
    }
}

impl StandardCombiner {
    fn fun_sim_avg(&self, m: &Matrix<f32>) -> f32 {
        let (rows, cols) = self.dim_f32(m);
        let row_maxes = self.row_maxes(m);
        let col_maxes = self.col_maxes(m);
        let mut nom = row_maxes.iter().sum::<f32>() / rows;
        nom += col_maxes.iter().sum::<f32>() / cols;

        nom / 2.0
    }

    fn fun_sim_max(&self, m: &Matrix<f32>) -> f32 {
        let (rows, cols) = self.dim_f32(m);
        let row_maxes = self.row_maxes(m);
        let col_maxes = self.col_maxes(m);

        (row_maxes.iter().sum::<f32>() / rows).max(col_maxes.iter().sum::<f32>() / cols)
    }

    fn bwa(&self, m: &Matrix<f32>) -> f32 {
        let (rows, cols) = self.dim_f32(m);
        let row_maxes = self.row_maxes(m);
        let col_maxes = self.col_maxes(m);

        (row_maxes.iter().sum::<f32>() + col_maxes.iter().sum::<f32>()) / (rows + cols)
    }
}

impl SimilarityCombiner for StandardCombiner {
    fn combine(&self, m: &Matrix<f32>) -> f32 {
        match self {
            StandardCombiner::FunSimAvg => self.fun_sim_avg(m),
            StandardCombiner::FunSimMax => self.fun_sim_max(m),
            StandardCombiner::Bwa => self.bwa(m),
        }
    }
}

/// calculate the Similarity score between two sets of HPO terms
pub struct GroupSimilarity<T, C> {
    combiner: C,
    similarity: T,
}

impl<T: Similarity, C: SimilarityCombiner> GroupSimilarity<T, C> {
    ///
    /// # Examples
    ///
    /// ```
    /// use hpo::similarity::GraphIc;
    /// use hpo::term::InformationContentKind;
    /// use hpo::similarity::{GroupSimilarity, StandardCombiner};
    ///
    /// // use Omim-based InformationContent for similarity calculation
    /// let graphic = GraphIc::new(InformationContentKind::Omim);
    ///
    /// // use the funSimAvg algorithm to combine the similarity scores
    /// let combiner = StandardCombiner::FunSimAvg;
    ///
    /// let sim = GroupSimilarity::new(combiner, graphic);
    /// ```
    ///
    pub fn new(combiner: C, similarity: T) -> Self {
        Self {
            combiner,
            similarity,
        }
    }

    /// calculates the similarity between two sets of terms
    pub fn calculate(&self, a: &HpoSet, b: &HpoSet) -> f32 {
        let mut v = Vec::with_capacity(a.len() * b.len());
        for t1 in a {
            for t2 in b {
                v.push(self.similarity.calculate(&t1, &t2));
            }
        }
        let m = Matrix::new(a.len(), b.len(), &v);
        self.combiner.calculate(&m)
    }
}

impl Default for GroupSimilarity<GraphIc, StandardCombiner> {
    fn default() -> Self {
        Self {
            combiner: StandardCombiner::default(),
            similarity: GraphIc::new(InformationContentKind::Omim),
        }
    }
}

/// This is a really weird way of converting a usize into a float but I
/// want to make sure the app crashes, so I don't want to use `as`.
fn usize_to_f32(n: usize) -> f32 {
    <usize as TryInto<u16>>::try_into(n)
        .expect("Matrix too large")
        .into()
}
