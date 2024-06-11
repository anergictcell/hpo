use tracing::debug;

use crate::annotations::GeneId;
use crate::stats::hypergeom::statrs::Hypergeometric;
use crate::stats::{f64_from_u64, Enrichment, SampleSet};
use crate::HpoTerm;

/// Calculates the hypergeometric enrichment of genes within the `set` compared to the `background`
///
/// # Examples
///
/// ```
/// use hpo::Ontology;
/// use hpo::{HpoSet, term::HpoGroup};
/// use hpo::stats::hypergeom::gene_enrichment;
///
/// let ontology = Ontology::from_binary("tests/example.hpo").unwrap();
///
/// let gene = ontology.gene_by_name("KRAS").unwrap();
/// let gene_hpo_set = gene.to_hpo_set(&ontology);
///
/// let mut enrichments = gene_enrichment(&ontology, &gene_hpo_set);
///
/// // the results are not sorted by default
/// enrichments.sort_by(|a, b| {
///         a.pvalue().partial_cmp(&b.pvalue()).unwrap()
/// });
/// assert!(enrichments.first().unwrap().pvalue() < enrichments.last().unwrap().pvalue());
/// assert!(enrichments.first().unwrap().enrichment() > enrichments.last().unwrap().enrichment());
/// ```
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
