//! Methods to calculate the Similarity between two terms or sets of terms
//!
//! Several methods and algorithms to calculate the similarity are already
//! provided in the library, but you can easily add your own as well.
//! The easiest way is to use the [`Builtins`] enum.
//!
//! # Examples
//!
//! ## Using built-in methods
//!
//! ```
//! use hpo::Ontology;
//! use hpo::similarity::{Builtins, Similarity};
//! use hpo::term::InformationContentKind;
//!
//! let ontology = Ontology::from_binary("tests/example.hpo").unwrap();
//! let term1 = ontology.hpo(12638u32.into()).unwrap();
//! let term2 = ontology.hpo(100547u32.into()).unwrap();
//!
//! let ic = Builtins::GraphIc(InformationContentKind::Omim);
//!
//! let similarity = ic.calculate(&term1, &term2);
//! println!("The termss {} and {} have a similarity of {}", term1.id(), term2.id(), similarity);
//! // ==> "The terms HP:0012638 and HP:0100547 have a similarity of 0.2704636"
//! ```
//!
//! ## Create a custom similarity algorithm
//! Creating you own similarity algorithm is as easy as implementing the
//! [Similarity](`crate::similarity::Similarity`) trait.
//!
//! ```
//! use hpo::{Ontology, HpoTerm};
//! use hpo::similarity::Similarity;
//!
//! struct Foo {}
//! impl Similarity for Foo {
//!     /// Calculate similarity based on length of the term names
//!     fn calculate(&self, a: &HpoTerm, b: &HpoTerm) -> f32 {
//!         return (a.name().len() - b.name().len()) as f32
//!     }
//! }
//!
//! let ontology = Ontology::from_binary("tests/example.hpo").unwrap();
//! let term1 = ontology.hpo(12638u32.into()).unwrap();
//! // ==> "Abnormal nervous system physiology"
//! let term2 = ontology.hpo(100547u32.into()).unwrap();
//! // ==> "Abnormal forebrain morphology"
//!
//! let ic = Foo{};
//!
//! let similarity = ic.calculate(&term1, &term2);
//! assert_eq!(similarity, 5.0);
//! ```

use std::cell::RefCell;
use std::collections::HashMap;

use crate::matrix::Matrix;
use crate::set::HpoSet;
use crate::term::InformationContentKind;
use crate::{HpoError, HpoResult, HpoTerm, HpoTermId};

pub mod defaults;
pub use defaults::{
    Distance, GraphIc, InformationCoefficient, Jc, Lin, Mutation, Relevance, Resnik,
};

/// Trait for similarity score calculation between 2 [`HpoTerm`]s
///
/// `hpo` already comes pre-loaded with several common and well established
/// similarity algorithms that implement the `Similarity` trait:
/// [Builtins](`crate::similarity::Builtins`)
///
/// ```
/// use hpo::{Ontology, HpoTerm};
/// use hpo::similarity::Similarity;
///
/// struct Foo {}
/// impl Similarity for Foo {
///     /// Calculate similarity based on length of the term names
///     fn calculate(&self, a: &HpoTerm, b: &HpoTerm) -> f32 {
///         return (a.name().len() - b.name().len()) as f32
///     }
/// }
///
/// let ontology = Ontology::from_binary("tests/example.hpo").unwrap();
/// let term1 = ontology.hpo(12638u32.into()).unwrap();
/// // ==> "Abnormal nervous system physiology"
/// let term2 = ontology.hpo(100547u32.into()).unwrap();
/// // ==> "Abnormal forebrain morphology"
///
/// let ic = Foo{};
///
/// let similarity = ic.calculate(&term1, &term2);
/// assert_eq!(similarity, 5.0);
/// ```
pub trait Similarity {
    /// calculates the actual similarity between term a and term b
    fn calculate(&self, a: &HpoTerm, b: &HpoTerm) -> f32;
}

/// This trait is needed to calculate the similarity between [`HpoSet`]s.
///
/// For similarity calculation between [`HpoSet`]s
/// the similarity scores must be combined to derive a single `f32` value
/// from a matrix of term - term similarities
///
/// `hpo` provides some default implementations of `SimilarityCombiner`:
/// [`StandardCombiner`](`crate::similarity::StandardCombiner`)
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
/// of 2 [`HpoSet`]s
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

impl TryFrom<&str> for StandardCombiner {
    type Error = HpoError;
    fn try_from(value: &str) -> Result<Self, Self::Error> {
        match value.to_lowercase().as_str() {
            "funsimavg" => Ok(StandardCombiner::FunSimAvg),
            "funsimmax" => Ok(StandardCombiner::FunSimMax),
            "bwa" => Ok(StandardCombiner::Bwa),
            _ => Err(HpoError::NotImplemented),
        }
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

/// calculate the Similarity score between two [`HpoSet`](`crate::HpoSet`)s
///
/// # Note
///
/// It is recommended to use the [`HpoSet::similarity`](`crate::HpoSet::similarity`)
/// method instead of creating a `GroupSimilarity` struct yourself.
///
/// # Examples
/// ## Using the preferred way
/// ```
/// use hpo::term::InformationContentKind;
/// use hpo::{Ontology, HpoSet};
/// use hpo::term::HpoGroup;
/// use hpo::similarity::{Builtins, StandardCombiner};
///
/// fn set1(ontology: &Ontology) -> HpoSet {
/// // ...
/// # let mut hpos = HpoGroup::new();
/// # hpos.insert(707u32.into());
/// # hpos.insert(12639u32.into());
/// # hpos.insert(12638u32.into());
/// # hpos.insert(818u32.into());
/// # hpos.insert(2715u32.into());
/// # HpoSet::new(ontology, hpos)
/// }
///
/// fn set2(ontology: &Ontology) -> HpoSet {
/// // ...
/// # let mut hpos = HpoGroup::new();
/// # hpos.insert(100547u32.into());
/// # hpos.insert(12638u32.into());
/// # hpos.insert(864u32.into());
/// # hpos.insert(25454u32.into());
/// # HpoSet::new(ontology, hpos)
/// }
///
/// let ontology = Ontology::from_binary("tests/example.hpo").unwrap();
/// let set_1 = set1(&ontology);
/// let set_2 = set2(&ontology);
///
/// let similarity = set_1.similarity(
///     &set_2,
///     Builtins::GraphIc(InformationContentKind::Omim),
///     StandardCombiner::default()
/// );
///
/// assert_eq!(similarity, 0.8177036);
/// ```
///
/// ## Using `GroupSimilarity` directly
///
/// ```
/// use hpo::term::InformationContentKind;
/// use hpo::{Ontology, HpoSet};
/// use hpo::term::HpoGroup;
/// use hpo::similarity::{Builtins, GroupSimilarity, StandardCombiner};
///
/// fn set1(ontology: &Ontology) -> HpoSet {
/// // ...
/// # let mut hpos = HpoGroup::new();
/// # hpos.insert(707u32.into());
/// # hpos.insert(12639u32.into());
/// # hpos.insert(12638u32.into());
/// # hpos.insert(818u32.into());
/// # hpos.insert(2715u32.into());
/// # HpoSet::new(ontology, hpos)
/// }
///
/// fn set2(ontology: &Ontology) -> HpoSet {
/// // ...
/// # let mut hpos = HpoGroup::new();
/// # hpos.insert(100547u32.into());
/// # hpos.insert(12638u32.into());
/// # hpos.insert(864u32.into());
/// # hpos.insert(25454u32.into());
/// # HpoSet::new(ontology, hpos)
/// }
///
/// let ontology = Ontology::from_binary("tests/example.hpo").unwrap();
/// let set_1 = set1(&ontology);
/// let set_2 = set2(&ontology);
///
///
/// let sim = GroupSimilarity::new(
///     StandardCombiner::FunSimAvg,
///     Builtins::GraphIc(InformationContentKind::Omim)
/// );
///
/// assert_eq!(sim.calculate(&set_1, &set_2), 0.8177036);
/// ```
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

/// Contains similarity methods for the standard built-in algorithms
///
/// For more details about each algorithm, check the [`defaults`] description.
///
/// # Examples
///
/// ```
/// use hpo::{Ontology, HpoSet};
/// use hpo::term::{InformationContentKind, HpoGroup};
/// use hpo::similarity::{Builtins, StandardCombiner};
///
/// fn set1(ontology: &Ontology) -> HpoSet {
/// // ...
/// # let mut hpos = HpoGroup::new();
/// # hpos.insert(707u32.into());
/// # hpos.insert(12639u32.into());
/// # hpos.insert(12638u32.into());
/// # hpos.insert(818u32.into());
/// # hpos.insert(2715u32.into());
/// # HpoSet::new(ontology, hpos)
/// }
///
/// fn set2(ontology: &Ontology) -> HpoSet {
/// // ...
/// # let mut hpos = HpoGroup::new();
/// # hpos.insert(100547u32.into());
/// # hpos.insert(12638u32.into());
/// # hpos.insert(864u32.into());
/// # hpos.insert(25454u32.into());
/// # HpoSet::new(ontology, hpos)
/// }
///
/// let ontology = Ontology::from_binary("tests/example.hpo").unwrap();
/// let set_1 = set1(&ontology);
/// let set_2 = set2(&ontology);
///
/// let sim_method = Builtins::GraphIc(InformationContentKind::Omim);
///
/// let similarity = set_1.similarity(
///     &set_2,
///     sim_method,
///     StandardCombiner::default()
/// );
/// ```
pub enum Builtins {
    /// [Distance](`Distance`) - based similarity
    Distance(InformationContentKind),
    /// [GraphIc](`GraphIc`) - based similarity
    GraphIc(InformationContentKind),
    /// [InformationCoefficient](`InformationCoefficient`) - based similarity
    InformationCoefficient(InformationContentKind),
    /// [Jiang & Conrath](`Jc`) - based similarity
    Jc(InformationContentKind),
    /// [Lin](`Lin`) - based similarity
    Lin(InformationContentKind),
    /// [Mutation](`Mutation`) - based similarity
    Mutation(InformationContentKind),
    /// [Relevance](`Relevance`) - based similarity
    Relevance(InformationContentKind),
    /// [Resnik](`Resnik`) - based similarity
    Resnik(InformationContentKind),
}

impl Builtins {
    /// Constructs a new `Builtins` struct from a `str`
    ///
    /// This method is useful to get a Similarity algorithm from a user provided string
    ///
    /// ```
    /// use hpo::term::InformationContentKind;
    /// use hpo::similarity::Builtins;
    ///
    /// let sim_method = Builtins::new("graphic", InformationContentKind::Omim);
    /// assert!(sim_method.is_ok());
    ///
    /// let sim_method = Builtins::new("does-not-exist", InformationContentKind::Omim);
    /// assert!(sim_method.is_err());
    /// ```
    ///
    /// # Errors
    ///
    /// Returns an [`HpoError::DoesNotExist`] error if no similary method with the given name exists
    pub fn new(method: &str, kind: InformationContentKind) -> HpoResult<Self> {
        match method.to_lowercase().as_str() {
            "graphic" => Ok(Self::GraphIc(kind)),
            "resnik" => Ok(Self::Resnik(kind)),
            "distance" | "dist" => Ok(Self::Distance(kind)),
            "informationcoefficient" | "ic" => Ok(Self::InformationCoefficient(kind)),
            "jc" | "jc2" => Ok(Self::Jc(kind)),
            "lin" => Ok(Self::Lin(kind)),
            "relevance" | "rel" => Ok(Self::Relevance(kind)),
            "mutation" | "mut" => Ok(Self::Mutation(kind)),
            _ => Err(HpoError::DoesNotExist),
        }
    }
}

impl Similarity for Builtins {
    fn calculate(&self, a: &HpoTerm, b: &HpoTerm) -> f32 {
        match self {
            Self::GraphIc(kind) => GraphIc::new(*kind).calculate(a, b),
            Self::Resnik(kind) => Resnik::new(*kind).calculate(a, b),
            Self::Distance(_) => Distance::new().calculate(a, b),
            Self::InformationCoefficient(kind) => {
                InformationCoefficient::new(*kind).calculate(a, b)
            }
            Self::Jc(kind) => Jc::new(*kind).calculate(a, b),
            Self::Lin(kind) => Lin::new(*kind).calculate(a, b),
            Self::Relevance(kind) => Relevance::new(*kind).calculate(a, b),
            Self::Mutation(kind) => Mutation::new(*kind).calculate(a, b),
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
