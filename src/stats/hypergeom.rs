//! Calculate the enrichment of a gene or disease in an `HpoSet` and the probability
//! of enrichment within the hypergeometric distribution.
//!
//! These methods are useful when you have clinical information of a patient
//! and want to see if the terms are enriched for some genes or diseases.
//!
//! # Examples
//!
//! ```no_run
//! use hpo::Ontology;
//! use hpo::{HpoSet, term::HpoGroup};
//! use hpo::stats::hypergeom::gene_enrichment;
//!
//! let ontology = Ontology::from_binary("tests/ontology.hpo").unwrap();
//!
//! let mut hpos = HpoGroup::new();
//! hpos.insert(2943u32.into());
//! hpos.insert(8458u32.into());
//! hpos.insert(100884u32.into());
//! hpos.insert(2944u32.into());
//! hpos.insert(2751u32.into());
//!
//! let patient_ci = HpoSet::new(&ontology, hpos);
//!
//! let mut enrichments = gene_enrichment(&ontology, &patient_ci);
//!     
//! // the results are not sorted by default
//! enrichments.sort_by(|a, b| {
//!         a.pvalue().partial_cmp(&b.pvalue()).unwrap()
//!     });
//!
//! for gene in enrichments {
//!     println!("{}\t{}\t({})", gene.id(), gene.pvalue(), gene.enrichment());
//! }
//! ```
//! You can also use it to find genes that have a similar phenotype to a gene.
//! Here we are creating an `HpoSet` from all [`HpoTerm`]s that are directly connected
//! to the gene `EZH2`. We're then checking which genes are enriched in the
//! directly and indirectly linked `HpoTerm`s.
//!
//! ```no_run
//! use hpo::Ontology;
//! use hpo::{HpoSet, term::HpoGroup};
//! use hpo::stats::hypergeom::gene_enrichment;
//!
//! let ontology = Ontology::from_binary("tests/ontology.hpo").unwrap();
//!
//! let gene = ontology.gene_by_name("EZH2").unwrap();
//! let gene_hpo_set = gene.to_hpo_set(&ontology);
//!
//! let mut enrichments = gene_enrichment(&ontology, &gene_hpo_set);
//!
//! // the results are not sorted by default
//! enrichments.sort_by(|a, b| {
//!         a.pvalue().partial_cmp(&b.pvalue()).unwrap()
//! });
//!
//! for gene in enrichments {
//!     println!("{}\t{}\t({})", gene.id(), gene.pvalue(), gene.enrichment());
//! }
//! ```

use log::debug;
use statrs::distribution::{DiscreteCDF, Hypergeometric};

use crate::annotations::{GeneId, OmimDiseaseId};
use crate::stats::{f64_from_u64, Enrichment, SampleSet};
use crate::HpoTerm;

/// Calculates the hypergeometric enrichment of genes within the `set` compared to the `background`
pub fn gene_enrichment<'a, T, U>(background: T, set: U) -> Vec<Enrichment<GeneId>>
where
    T: IntoIterator<Item = HpoTerm<'a>>,
    U: IntoIterator<Item = HpoTerm<'a>>,
{
    fn inner_gene_enrichment(
        background: &SampleSet<GeneId>,
        sample_set: &SampleSet<GeneId>,
    ) -> Vec<Enrichment<GeneId>> {
        let mut res = Vec::new();
        for (gene, observed_successes) in sample_set {
            if observed_successes == 0 {
                debug!("Skipping {}", gene);
                continue;
            }
            let successes = background
                .get(&gene)
                .expect("gene must be present in background set");
            let hyper = Hypergeometric::new(
                // Total number of HPOTerms in the Ontology
                // ==> population
                background.len(),
                // Number of terms in the Ontology that are associated to the gene
                // ==> successes
                *successes,
                // Number of terms in the set
                // ==> draws
                sample_set.len(),
            )
            .expect("the set must not be larger than the ontology");

            // subtracting 1, because we want to test including observed_successes
            // e.g. "7 or more", but sf by default calculates "more than 7"
            let pvalue = hyper.sf(observed_successes - 1);
            let enrichment = (f64_from_u64(observed_successes) / f64_from_u64(sample_set.len()))
                / (f64_from_u64(*successes) / f64_from_u64(background.len()));

            res.push(Enrichment::gene(
                gene,
                pvalue,
                observed_successes,
                enrichment,
            ));
            debug!(
                "Gene:{}\tPopulation: {}, Successes: {}, Draws: {}, Observed: {}",
                gene,
                background.len(),
                successes,
                sample_set.len(),
                observed_successes
            );
        }
        res
    }

    let background = SampleSet::gene(background);
    let sample_set = SampleSet::gene(set);
    inner_gene_enrichment(&background, &sample_set)
}

/// Calculates the hypergeometric enrichment of diseases within the `set` compared to the `background`
pub fn disease_enrichment<'a, T, U>(background: T, set: U) -> Vec<Enrichment<OmimDiseaseId>>
where
    T: IntoIterator<Item = HpoTerm<'a>>,
    U: IntoIterator<Item = HpoTerm<'a>>,
{
    fn inner_disease_enrichment(
        background: &SampleSet<OmimDiseaseId>,
        sample_set: &SampleSet<OmimDiseaseId>,
    ) -> Vec<Enrichment<OmimDiseaseId>> {
        let mut res = Vec::new();
        for (disease, observed_successes) in sample_set {
            if observed_successes == 0 {
                debug!("Skipping {}", disease);
                continue;
            }
            let successes = background
                .get(&disease)
                .expect("disease must be present in background set");
            let hyper = Hypergeometric::new(
                // Total number of HPOTerms in the Ontology
                // ==> population
                background.len(),
                // Number of terms in the Ontology that are associated to the disease
                // ==> successes
                *successes,
                // Number of terms in the set
                // ==> draws
                sample_set.len(),
            )
            .expect("the set must not be larger than the background");
            // subtracting 1, because we want to test including observed_successes
            // e.g. "7 or more", but sf by default calculates "more than 7"
            let pvalue = hyper.sf(observed_successes - 1);
            let enrichment = (f64_from_u64(observed_successes) / f64_from_u64(sample_set.len()))
                / (f64_from_u64(*successes) / f64_from_u64(background.len()));
            res.push(Enrichment::disease(
                disease,
                pvalue,
                observed_successes,
                enrichment,
            ));
            debug!(
                "Disease:{}\tPopulation: {}, Successes: {}, Draws: {}, Observed: {}",
                disease,
                background.len(),
                successes,
                sample_set.len(),
                observed_successes
            );
        }
        res
    }

    let background = SampleSet::disease(background);
    let sample_set = SampleSet::disease(set);
    inner_disease_enrichment(&background, &sample_set)
}
