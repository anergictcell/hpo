use tracing::debug;

use crate::annotations::OmimDiseaseId;
use crate::stats::hypergeom::statrs::Hypergeometric;
use crate::stats::{f64_from_u64, Enrichment, SampleSet};
use crate::HpoTerm;

/// Calculates the hypergeometric enrichment of diseases within the `set` compared to the `background`
///
/// # Examples
///
/// ```
/// use hpo::Ontology;
/// use hpo::{HpoSet, term::HpoGroup};
/// use hpo::stats::hypergeom::disease_enrichment;
///
/// let ontology = Ontology::from_binary("tests/example.hpo").unwrap();
///
/// let gene = ontology.gene_by_name("KRAS").unwrap();
/// let gene_hpo_set = gene.to_hpo_set(&ontology);
///
/// let mut enrichments = disease_enrichment(&ontology, &gene_hpo_set);
///
/// // the results are not sorted by default
/// enrichments.sort_by(|a, b| {
///         a.pvalue().partial_cmp(&b.pvalue()).unwrap()
/// });
/// assert!(enrichments.first().unwrap().pvalue() < enrichments.last().unwrap().pvalue());
/// assert!(enrichments.first().unwrap().enrichment() > enrichments.last().unwrap().enrichment());
/// ```
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
            res.push(Enrichment::annotation(
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
