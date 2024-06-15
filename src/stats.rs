//! Statistical analyses for [`HpoTerm`] annotations and enrichment
//!
//! This module contains methods to calculate the enrichment of Genes or Diseases
//! within sets of `HpoTerm`s, such as a patient's clinical information.
//!
//! It can also be used to identify genes or diseases with a similar phenotype,
//! even comparing genes to identify novel gene-interaction-networks.
//!
//! At the moment, `hpo` provides only the hypergeometric enrichment analysis using
//! the survival function.
//!
//! In addition, it provides methods for [hierarchical clustering](`linkage::Linkage`) of `HpoSet`s.

use std::collections::HashMap;
use std::marker::PhantomData;

use crate::annotations::{AnnotationId, Disease, GeneId, OmimDiseaseId, OrphaDiseaseId};
use crate::HpoTerm;

pub mod hypergeom;
mod linkage;
pub use linkage::cluster;
pub use linkage::Linkage;

/// The fold enrichment and p-value for an enriched Gene or Disease
///
/// [`Enrichment`] is returned from statistics enrichment methods, such as
/// [`hypergeom::gene_enrichment`], [`hypergeom::omim_disease_enrichment`] and [`hypergeom::orpha_disease_enrichment`].
///
#[derive(Debug)]
pub struct Enrichment<T> {
    annotation: T,
    pvalue: f64,
    count: u64,
    enrichment: f64,
}

impl<T> Enrichment<T> {
    /// Returns the p-value of the enrichment
    ///
    /// The p-value indicates the probability that the enrichment
    /// occured by chance
    pub fn pvalue(&self) -> f64 {
        self.pvalue
    }

    /// Returns the fold enrichment over the background population
    pub fn enrichment(&self) -> f64 {
        self.enrichment
    }

    /// Returns the ID of the enrichment item, most likely `GeneId` or `OmimDiseaseId`
    pub fn id(&self) -> &T {
        &self.annotation
    }

    /// Returns the number of items in the enrichment set
    pub fn count(&self) -> u64 {
        self.count
    }
}

impl Enrichment<GeneId> {
    /// Constructs an `Enrichment` for `GeneId`s
    pub fn gene(gene: GeneId, pvalue: f64, count: u64, enrichment: f64) -> Self {
        Self {
            annotation: gene,
            pvalue,
            count,
            enrichment,
        }
    }
}

impl Enrichment<OmimDiseaseId> {
    /// Constructs an `Enrichment` for `OmimDiseaseId`s
    pub fn disease(disease: OmimDiseaseId, pvalue: f64, count: u64, enrichment: f64) -> Self {
        Self {
            annotation: disease,
            pvalue,
            count,
            enrichment,
        }
    }
}

impl<T: AnnotationId> Enrichment<T> {
    /// Constructs an `Enrichment` for any annotated item (Gene, Disease)
    pub fn annotation(annotation: T, pvalue: f64, count: u64, enrichment: f64) -> Self {
        Self {
            annotation,
            pvalue,
            count,
            enrichment,
        }
    }
}

struct SampleSet<T> {
    /// The total number of `HpoTerms` in the full set
    size: u64,
    /// A map containing the counts of each Gene/Disease in the `SampleSet`
    counts: HashMap<u32, u64>,
    phantom: PhantomData<T>,
}

fn calculate_counts<
    'a,
    U: FnMut(HpoTerm<'a>) -> IT,
    I: IntoIterator<Item = HpoTerm<'a>>,
    IT: IntoIterator<Item = u32>,
>(
    terms: I,
    mut iter: U,
) -> (u64, HashMap<u32, u64>) {
    let mut size = 0u64;
    let mut counts: HashMap<u32, u64> = HashMap::new();
    for term in terms {
        size += 1;
        for id in iter(term) {
            counts
                .entry(id)
                .and_modify(|count| *count += 1)
                .or_insert(1);
        }
    }
    (size, counts)
}

impl<'a> SampleSet<GeneId> {
    /// Constructs a new [`SampleSet`] with gene counts from an iterator of [`HpoTerm`]s
    pub fn gene<I: IntoIterator<Item = HpoTerm<'a>>>(terms: I) -> Self {
        let term2geneid = |term: HpoTerm<'a>| term.genes().map(|d| d.id().as_u32());

        let (size, counts) = calculate_counts(terms, term2geneid);
        Self {
            size,
            counts,
            phantom: PhantomData,
        }
    }
}

impl<'a> SampleSet<OmimDiseaseId> {
    /// Constructs a new `SampleSet` with disease counts from an iterator of [`HpoTerm`]s
    pub fn omim_disease<I: IntoIterator<Item = HpoTerm<'a>>>(terms: I) -> Self {
        let term2omimid = |term: HpoTerm<'a>| term.omim_diseases().map(|d| d.id().as_u32());
        let (size, counts) = calculate_counts(terms, term2omimid);
        Self {
            size,
            counts,
            phantom: PhantomData,
        }
    }
}

impl<'a> SampleSet<OrphaDiseaseId> {
    /// Constructs a new `SampleSet` with disease counts from an iterator of [`HpoTerm`]s
    pub fn orpha_disease<I: IntoIterator<Item = HpoTerm<'a>>>(terms: I) -> Self {
        let term2omimid = |term: HpoTerm<'a>| term.orpha_diseases().map(|d| d.id().as_u32());
        let (size, counts) = calculate_counts(terms, term2omimid);
        Self {
            size,
            counts,
            phantom: PhantomData,
        }
    }
}

impl<T: AnnotationId> SampleSet<T> {
    /// Returns the total number of [`HpoTerm`]s in the [`SampleSet`]
    fn len(&self) -> u64 {
        self.size
    }

    /// Returns whether the set is empty
    ///
    /// Returns `true` if the set is empty, `false` otherwise
    #[allow(dead_code)]
    fn is_empty(&self) -> bool {
        self.size == 0
    }
    /// The number of terms in the `SampleSet` that are connected to the given gene or disease
    ///
    /// Terms can be directly or indirectly connected to the gene/disease.
    ///
    /// Returns `None` if the key is not present
    fn get(&self, key: &T) -> Option<&u64> {
        self.counts.get(&key.as_u32())
    }

    /// Returns the frequency of terms that are connected to the given gene or disease
    ///
    /// It divides the number of terms that are connected the gene/disease
    /// by the total number of terms in the `SampleSet`.
    ///
    /// Returns `None` if the key is not present
    #[allow(dead_code)]
    fn frequency(&self, key: &T) -> Option<f64> {
        self.counts
            .get(&key.as_u32())
            .map(|count| f64_from_u64(*count) / f64_from_u64(self.size))
    }

    /// An iterator of [`SampleSet::frequency`] values, along with their key
    #[allow(dead_code)]
    fn frequencies(&self) -> Frequencies<T> {
        Frequencies::new(self.counts.iter(), self.size, self.phantom)
    }
}

impl<'a, T: AnnotationId> IntoIterator for &'a SampleSet<T> {
    type Item = (T, u64);
    type IntoIter = Counts<'a, T>;
    fn into_iter(self) -> Self::IntoIter {
        Counts::new(self.counts.iter(), self.phantom)
    }
}

/// An iterator of [`SampleSet::frequency`] values, along with their key
struct Frequencies<'a, K> {
    inner: std::collections::hash_map::Iter<'a, u32, u64>,
    total: u64,
    phantom: PhantomData<K>,
}

impl<'a, K> Frequencies<'a, K> {
    pub fn new(
        inner: std::collections::hash_map::Iter<'a, u32, u64>,
        total: u64,
        phantom: PhantomData<K>,
    ) -> Self {
        Self {
            inner,
            total,
            phantom,
        }
    }
}

impl<K: AnnotationId> Iterator for Frequencies<'_, K> {
    type Item = (K, f64);
    fn next(&mut self) -> Option<Self::Item> {
        self.inner
            .next()
            .map(|(k, v)| (K::from(*k), f64_from_u64(*v) / f64_from_u64(self.total)))
    }
}

/// An iterator of [`SampleSet::frequency`] values, along with their key
struct Counts<'a, K> {
    inner: std::collections::hash_map::Iter<'a, u32, u64>,
    phantom: PhantomData<K>,
}

impl<'a, K> Counts<'a, K> {
    pub fn new(
        inner: std::collections::hash_map::Iter<'a, u32, u64>,
        phantom: PhantomData<K>,
    ) -> Self {
        Self { inner, phantom }
    }
}

impl<K: AnnotationId> Iterator for Counts<'_, K> {
    type Item = (K, u64);
    fn next(&mut self) -> Option<Self::Item> {
        self.inner.next().map(|(k, v)| (K::from(*k), *v))
    }
}

/// We have to frequently do divisions starting with u64 values
/// and need to return f64 values. To ensure some kind of safety
/// we use this method to panic in case of overflows.
fn f64_from_u64(n: u64) -> f64 {
    let intermediate: u32 = n
        .try_into()
        .expect("cannot safely create f64 from large u64");
    intermediate.into()
}

/// We have to frequently do divisions starting with u64 values
/// and need to return f64 values. To ensure some kind of safety
/// we use this method to panic in case of overflows.
fn f64_from_usize(n: usize) -> f64 {
    let intermediate: u32 = n
        .try_into()
        .expect("cannot safely create f64 from large u64");
    intermediate.into()
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn iterate_frequencies() {
        let mut map = HashMap::new();
        map.insert(12u32, 12u64);
        map.insert(21u32, 21u64);

        let mut iter: Frequencies<'_, OmimDiseaseId> = Frequencies::new(map.iter(), 3, PhantomData);
        match iter.next() {
            Some((key, x)) if key == OmimDiseaseId::from(12) => {
                assert!((x - 4.0).abs() < f64::EPSILON);
            }
            Some((key, x)) if key == OmimDiseaseId::from(21) => {
                assert!((x - 7.0).abs() < f64::EPSILON);
            }
            _ => panic!("invalid"),
        }
        match iter.next() {
            Some((key, x)) if key == OmimDiseaseId::from(12) => {
                assert!((x - 4.0).abs() < f64::EPSILON);
            }
            Some((key, x)) if key == OmimDiseaseId::from(21) => {
                assert!((x - 7.0).abs() < f64::EPSILON);
            }
            _ => panic!("invalid"),
        }
        assert!(iter.next().is_none());
    }

    #[test]
    fn iterate_counts() {
        let mut map = HashMap::new();
        map.insert(12u32, 12u64);
        map.insert(21u32, 21u64);

        let mut iter: Counts<'_, OmimDiseaseId> = Counts::new(map.iter(), PhantomData);
        match iter.next() {
            Some((key, x)) if key == OmimDiseaseId::from(12) => assert_eq!(x, 12),
            Some((key, x)) if key == OmimDiseaseId::from(21) => assert_eq!(x, 21),
            _ => panic!("invalid"),
        }
        match iter.next() {
            Some((key, x)) if key == OmimDiseaseId::from(12) => assert_eq!(x, 12),
            Some((key, x)) if key == OmimDiseaseId::from(21) => assert_eq!(x, 21),
            _ => panic!("invalid"),
        }
        assert!(iter.next().is_none());
    }
}
