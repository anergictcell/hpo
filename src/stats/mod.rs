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

use std::collections::HashMap;
use std::hash::Hash;

use crate::annotations::{GeneId, OmimDiseaseId};
use crate::HpoTerm;

pub mod hypergeom;

/// The fold enrichment and p-value for an enriched Gene or Disease
///
/// [`Enrichment`] is returned from statistics enrichment methods, such as
/// [`hypergeom::gene_enrichment`] and [`hypergeom::disease_enrichment`]
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

struct SampleSet<T> {
    /// The total number of HpoTerms in the full set
    size: u64,
    /// A map containing the counts of each Gene/Disease in the SampleSet
    counts: HashMap<T, u64>,
}

impl<'a> SampleSet<GeneId> {
    /// Constructs a new [`SampleSet`] with gene counts from an iterator of [`HpoTerm`]s
    pub fn gene<I: IntoIterator<Item = HpoTerm<'a>>>(iter: I) -> Self {
        let mut size = 0u64;
        let mut counts: HashMap<GeneId, u64> = HashMap::new();
        for term in iter {
            size += 1;
            for gene in term.genes() {
                counts
                    .entry(*gene.id())
                    .and_modify(|count| *count += 1)
                    .or_insert(1);
            }
        }
        Self { size, counts }
    }
}

impl<'a> SampleSet<OmimDiseaseId> {
    /// Constructs a new `SampleSet` with disease counts from an iterator of [`HpoTerm`]s
    pub fn disease<I: IntoIterator<Item = HpoTerm<'a>>>(ontology: I) -> Self {
        let mut size = 0u64;
        let mut counts: HashMap<OmimDiseaseId, u64> = HashMap::new();
        for term in ontology {
            size += 1;
            for disease in term.omim_diseases() {
                counts
                    .entry(*disease.id())
                    .and_modify(|count| *count += 1)
                    .or_insert(1);
            }
        }
        Self { size, counts }
    }
}

impl<T> SampleSet<T> {
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
}

impl<T: Hash + std::cmp::Eq> SampleSet<T> {
    /// The number of terms in the `SampleSet` that are connected to the given gene or disease
    ///
    /// Terms can be directly or indirectly connected to the gene/disease.
    ///
    /// Returns `None` if the key is not present
    fn get(&self, key: &T) -> Option<&u64> {
        self.counts.get(key)
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
            .get(key)
            .map(|count| f64_from_u64(*count) / f64_from_u64(self.size))
    }

    /// An iterator of [`SampleSet::frequency`] values, along with their key
    #[allow(dead_code)]
    fn frequencies(&self) -> Frequencies<T> {
        Frequencies::new(self.counts.iter(), self.size)
    }
}

impl<'a, T: Clone> IntoIterator for &'a SampleSet<T> {
    type Item = (T, u64);
    type IntoIter = Counts<'a, T>;
    fn into_iter(self) -> Self::IntoIter {
        Counts::new(self.counts.iter())
    }
}

/// An iterator of [`SampleSet::frequency`] values, along with their key
struct Frequencies<'a, K> {
    inner: std::collections::hash_map::Iter<'a, K, u64>,
    total: u64,
}

impl<'a, K> Frequencies<'a, K> {
    pub fn new(inner: std::collections::hash_map::Iter<'a, K, u64>, total: u64) -> Self {
        Self { inner, total }
    }
}

impl<K: Clone> Iterator for Frequencies<'_, K> {
    type Item = (K, f64);
    fn next(&mut self) -> Option<Self::Item> {
        self.inner
            .next()
            .map(|(k, v)| (k.clone(), f64_from_u64(*v) / f64_from_u64(self.total)))
    }
}

/// An iterator of [`SampleSet::frequency`] values, along with their key
struct Counts<'a, K> {
    inner: std::collections::hash_map::Iter<'a, K, u64>,
}

impl<'a, K> Counts<'a, K> {
    pub fn new(inner: std::collections::hash_map::Iter<'a, K, u64>) -> Self {
        Self { inner }
    }
}

impl<K: Clone> Iterator for Counts<'_, K> {
    type Item = (K, u64);
    fn next(&mut self) -> Option<Self::Item> {
        self.inner.next().map(|(k, v)| (k.clone(), *v))
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn iterate_frequencies() {
        let mut map = HashMap::new();
        map.insert(String::from("foo"), 12u64);
        map.insert(String::from("bar"), 21u64);

        let mut iter = Frequencies::new(map.iter(), 3);
        match iter.next() {
            Some((key, x)) if key == "foo" => assert_eq!(x, 4.0),
            Some((key, x)) if key == "bar" => assert_eq!(x, 7.0),
            _ => panic!("invalid"),
        }
        match iter.next() {
            Some((key, x)) if key == "foo" => assert_eq!(x, 4.0),
            Some((key, x)) if key == "bar" => assert_eq!(x, 7.0),
            _ => panic!("invalid"),
        }
        assert!(iter.next().is_none());
    }

    #[test]
    fn iterate_counts() {
        let mut map = HashMap::new();
        map.insert(String::from("foo"), 12u64);
        map.insert(String::from("bar"), 21u64);

        let mut iter = Counts::new(map.iter());
        match iter.next() {
            Some((key, x)) if key == "foo" => assert_eq!(x, 12),
            Some((key, x)) if key == "bar" => assert_eq!(x, 21),
            _ => panic!("invalid"),
        }
        match iter.next() {
            Some((key, x)) if key == "foo" => assert_eq!(x, 12),
            Some((key, x)) if key == "bar" => assert_eq!(x, 21),
            _ => panic!("invalid"),
        }
        assert!(iter.next().is_none());
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
