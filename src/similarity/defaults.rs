//! Contains implementation for some Similarity algorithms that are
//! considered to be default implementations.
//!
//! All of the algorithms can also be accessed via [`crate::similarity::Builtins`]

use crate::similarity::{usize_to_f32, Similarity};
use crate::term::InformationContentKind;
use crate::HpoTerm;

/// Graph based Information coefficient similarity
///
/// For a detailed description see [Deng Y, et. al., PLoS One, (2015)](https://pubmed.ncbi.nlm.nih.gov/25664462/)
pub struct GraphIc {
    kind: InformationContentKind,
}

impl GraphIc {
    /// Constructs a new struct to calculate GraphIC based similarity scores
    /// between two terms
    ///
    /// # Examples
    ///
    /// ```
    /// use hpo::similarity::GraphIc;
    /// use hpo::term::InformationContentKind;
    ///
    /// // use Omim-based InformationContent for similarity calculation
    /// let graphic = GraphIc::new(InformationContentKind::Omim);
    /// ```
    ///
    pub fn new(kind: InformationContentKind) -> Self {
        Self { kind }
    }
}

impl Similarity for GraphIc {
    fn calculate(&self, a: &HpoTerm, b: &HpoTerm) -> f32 {
        if a.id() == b.id() {
            return 1.0;
        }

        let ic_union: f32 = a
            .all_union_ancestors(b)
            .map(|p| p.information_content().get_kind(&self.kind))
            .sum();

        if ic_union == 0.0 {
            return 0.0;
        }

        let ic_common: f32 = a
            .all_common_ancestors(b)
            .map(|p| p.information_content().get_kind(&self.kind))
            .sum();

        ic_common / ic_union
    }
}

/// Similarity score from Resnik
///
/// For a detailed description see [Resnik P, Proceedings of the 14th IJCAI, (1995)](https://www.ijcai.org/Proceedings/95-1/Papers/059.pdf)
pub struct Resnik {
    kind: InformationContentKind,
}

impl Resnik {
    /// Constructs a new struct to calculate the Resnik based similarity scores
    /// between two terms
    ///
    /// # Examples
    ///
    /// ```
    /// use hpo::similarity::Resnik;
    /// use hpo::term::InformationContentKind;
    ///
    /// // use Omim-based InformationContent for similarity calculation
    /// let resnik = Resnik::new(InformationContentKind::Omim);
    /// ```
    ///
    pub fn new(kind: InformationContentKind) -> Self {
        Self { kind }
    }
}

impl Similarity for Resnik {
    fn calculate(&self, a: &HpoTerm, b: &HpoTerm) -> f32 {
        a.all_common_ancestors(b)
            .map(|term| term.information_content().get_kind(&self.kind))
            .fold(0.0, |max, term| if term > max { term } else { max })
    }
}

/// Similarity score from Lin
///
/// For a detailed description see [Lin D, Proceedings of the 15th ICML, (1998)](https://dl.acm.org/doi/10.5555/645527.657297)
pub struct Lin {
    kind: InformationContentKind,
}

impl Lin {
    /// Constructs a new struct to calculate the Lin based similarity scores
    /// between two terms
    ///
    /// # Examples
    ///
    /// ```
    /// use hpo::similarity::Lin;
    /// use hpo::term::InformationContentKind;
    ///
    /// // use Omim-based InformationContent for similarity calculation
    /// let lin = Lin::new(InformationContentKind::Omim);
    /// ```
    ///
    pub fn new(kind: InformationContentKind) -> Self {
        Self { kind }
    }
}

impl Similarity for Lin {
    fn calculate(&self, a: &HpoTerm, b: &HpoTerm) -> f32 {
        let ic_combined = a.information_content().get_kind(&self.kind)
            + b.information_content().get_kind(&self.kind);

        if ic_combined == 0.0 {
            return 0.0;
        }

        let resnik = Resnik::new(self.kind).calculate(a, b);

        2.0 * resnik / ic_combined
    }
}

/// Similarity score from Jiang & Conrath
///
/// For a detailed description see [Jiang J, Conrath D, ROCLING X, (1997)](https://aclanthology.org/O97-1002.pdf)
///
/// # Note
///
/// This algorithm is an implementation as described in the paper cited above. It is different
/// from the `JC` implementation in the `HPOSim` R library. It is identical to the `JC2`
/// implementation in [`PyHPO`](https://pypi.org/project/pyhpo/)
pub struct Jc {
    kind: InformationContentKind,
}

impl Jc {
    /// Constructs a new struct to calculate the Jiang & Conrath based similarity scores
    /// between two terms
    ///
    /// # Examples
    ///
    /// ```
    /// use hpo::similarity::Jc;
    /// use hpo::term::InformationContentKind;
    ///
    /// // use Omim-based InformationContent for similarity calculation
    /// let jc = Jc::new(InformationContentKind::Omim);
    /// ```
    ///
    pub fn new(kind: InformationContentKind) -> Self {
        Self { kind }
    }
}

impl Similarity for Jc {
    fn calculate(&self, a: &HpoTerm, b: &HpoTerm) -> f32 {
        if a.id() == b.id() {
            return 1.0;
        }

        let ic_combined = a.information_content().get_kind(&self.kind)
            + b.information_content().get_kind(&self.kind);

        let resnik = Resnik::new(self.kind).calculate(a, b);

        1.0 - (ic_combined - 2.0 * resnik)
    }
}

/// Relevance Similarity score from Schlicker
///
/// For a detailed description see [Schlicker A, et.al., BMC Bioinformatics, (2006)](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-7-302)
///
pub struct Relevance {
    kind: InformationContentKind,
}

impl Relevance {
    /// Constructs a new struct to calculate the Schlicker based similarity scores
    /// between two terms
    ///
    /// # Examples
    ///
    /// ```
    /// use hpo::similarity::Relevance;
    /// use hpo::term::InformationContentKind;
    ///
    /// // use Omim-based InformationContent for similarity calculation
    /// let rel = Relevance::new(InformationContentKind::Omim);
    /// ```
    ///
    pub fn new(kind: InformationContentKind) -> Self {
        Self { kind }
    }
}

impl Similarity for Relevance {
    fn calculate(&self, a: &HpoTerm, b: &HpoTerm) -> f32 {
        let resnik = Resnik::new(self.kind).calculate(a, b);
        let lin = Lin::new(self.kind).calculate(a, b);

        lin * (1.0 - (resnik * -1.0).exp())
    }
}

/// Information Coefficient Similarity score from Li
///
/// For a detailed description see [Li B, et. al., arXiv, (2010)](https://arxiv.org/abs/1001.0958)
///
pub struct InformationCoefficient {
    kind: InformationContentKind,
}

impl InformationCoefficient {
    /// Constructs a new struct to calculate the Jiang & Conrath based similarity scores
    /// between two terms
    ///
    /// # Examples
    ///
    /// ```
    /// use hpo::similarity::InformationCoefficient;
    /// use hpo::term::InformationContentKind;
    ///
    /// // use Omim-based InformationContent for similarity calculation
    /// let ic = InformationCoefficient::new(InformationContentKind::Omim);
    /// ```
    ///
    pub fn new(kind: InformationContentKind) -> Self {
        Self { kind }
    }
}

impl Similarity for InformationCoefficient {
    fn calculate(&self, a: &HpoTerm, b: &HpoTerm) -> f32 {
        let resnik = Resnik::new(self.kind).calculate(a, b);
        let lin = Lin::new(self.kind).calculate(a, b);

        lin * (1.0 - (1.0 / (1.0 + resnik)))
    }
}

/// Similarity score based on distance between terms
#[derive(Default)]
pub struct Distance {}

impl Distance {
    /// Constructs a new struct to calculate the distance based similarity scores
    /// between two terms
    ///
    /// # Examples
    ///
    /// ```
    /// use hpo::similarity::Distance;
    /// use hpo::term::InformationContentKind;
    ///
    /// let dist = Distance::new();
    /// ```
    ///
    pub fn new() -> Self {
        Self::default()
    }
}

impl Similarity for Distance {
    fn calculate(&self, a: &HpoTerm, b: &HpoTerm) -> f32 {
        a.distance_to_term(b)
            .map_or(0.0, |n| 1.0 / (usize_to_f32(n) + 1.0))
    }
}
