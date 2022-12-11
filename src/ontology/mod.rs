//! Ontology is the heart and main interface of the full crate.
//!
//! The [`Ontology`] struct holds all information about the ontology
//! and the ownership of all [`HpoTerm`]s, [`Gene`]s and [`OmimDisease`]s.

use std::collections::HashMap;
use std::ops::BitOr;
use std::path::Path;

use crate::annotations::{Gene, GeneId};
use crate::annotations::{OmimDisease, OmimDiseaseId};
use crate::term::internal::HpoTermInternal;
use crate::term::HpoTerm;
use crate::HpoTermId;
use crate::OntologyResult;
use crate::{parser, HpoParents};

use core::fmt::Debug;

mod termarena;
use termarena::Arena;

/// Main API interface that owns all data
///
/// It is recommended to use the public methods to build the ontology
/// from standard annotation data from Jax. You will need to download
/// the data from [HPO](https://hpo.jax.org/) itself.
#[derive(Default)]
pub struct Ontology {
    hpo_terms: Arena,
    genes: HashMap<GeneId, Gene>,
    omim_diseases: HashMap<OmimDiseaseId, OmimDisease>,
}

impl Debug for Ontology {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Ontology with {} terns", self.hpo_terms.len())
    }
}

/// Public API of the Ontology
///
/// Those methods are all safe to use
impl Ontology {
    /// Initialize the [`Ontology`] from data provided by [Jax HPO](https://hpo.jax.org/)
    ///
    /// You must download:
    ///
    /// - Actual OBO data: [hp.obo](https://hpo.jax.org/app/data/ontology)
    /// - Links between HPO and OMIM diseases: [phenotype.hpoa](https://hpo.jax.org/app/data/annotations)
    /// - Links between HPO and Genes: [phenotype_to_genes.txt](http://purl.obolibrary.org/obo/hp/hpoa/phenotype_to_genes.txt)
    ///
    /// and then specify the folder where the data is stored.
    ///
    /// # Note
    ///
    /// It is quite likely that the method signature will change and instead
    /// return a `Result`
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use hpo::Ontology;
    /// use hpo::HpoTermId;
    ///
    /// let ontology = Ontology::from_standard("./example_data/");
    ///
    /// assert!(ontology.len() > 15_000);
    ///
    /// let absent_term = HpoTermId::try_from("HP:9999999").unwrap();
    /// assert!(ontology.hpo(&absent_term).is_none());
    ///
    /// let present_term = HpoTermId::try_from("HP:0000001").unwrap();
    /// let root_term = ontology.hpo(&present_term).unwrap();
    /// assert_eq!(root_term.name(), "Phenotypical abnormality");
    /// ```
    ///
    pub fn from_standard(folder: &str) -> Self {
        let mut ont = Ontology::default();
        let path = Path::new(folder);
        let obo = path.join(crate::OBO_FILENAME);
        let gene = path.join(crate::GENE_FILENAME);
        let disease = path.join(crate::DISEASE_FILENAME);
        parser::load_from_standard_files(&obo, &gene, &disease, &mut ont);
        ont
    }

    /// Returns the number of HPO-Terms in the Ontology
    pub fn len(&self) -> usize {
        self.hpo_terms.len()
    }

    /// Returns `true` if the Ontology holds no HPO-Terms
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Returns the [`HpoTerm`] of the provided [`HpoTermId`]
    ///
    /// If no such term is present in the Ontolgy, `None` is returned
    pub fn hpo(&self, term_id: &HpoTermId) -> Option<HpoTerm> {
        HpoTerm::try_new(self, term_id).ok()
    }

    /// Returns an Iterator of all [`HpoTerm`]s from the Ontology
    pub fn hpos(&self) -> OntologyIterator {
        OntologyIterator {
            inner: self.hpo_terms.values().iter(),
            ontology: self,
        }
    }

    /// Return a reference to the [`Gene`] of the provided [`GeneId`]
    ///
    /// If no such gene is present, `None` is returned
    pub fn gene(&self, gene_id: &GeneId) -> Option<&Gene> {
        self.genes.get(gene_id)
    }

    /// Return a mutable reference to the [`Gene`] of the provided [`GeneId`]
    ///
    /// If no such gene is present, `None` is returned
    pub fn gene_mut(&mut self, gene_id: &GeneId) -> Option<&mut Gene> {
        self.genes.get_mut(gene_id)
    }

    /// Returns an Iterator of all [`Gene`]s from the Ontology
    pub fn genes(&self) -> std::collections::hash_map::Values<'_, GeneId, Gene> {
        self.genes.values()
    }

    /// Return a reference to the [`OmimDisease`] of the provided [`OmimDiseaseId`]
    ///
    /// If no such disease is present, `None` is returned
    pub fn omim_disease(&self, omim_disease_id: &OmimDiseaseId) -> Option<&OmimDisease> {
        self.omim_diseases.get(omim_disease_id)
    }

    /// Return a mutable reference to the [`OmimDisease`] of the provided [`OmimDiseaseId`]
    ///
    /// If no such disease is present, `None` is returned
    pub fn omim_disease_mut(
        &mut self,
        omim_disease_id: &OmimDiseaseId,
    ) -> Option<&mut OmimDisease> {
        self.omim_diseases.get_mut(omim_disease_id)
    }

    /// Returns an Iterator of all [`OmimDisease`]s from the Ontology
    pub fn omim_diseases(
        &self,
    ) -> std::collections::hash_map::Values<'_, OmimDiseaseId, OmimDisease> {
        self.omim_diseases.values()
    }

    pub(crate) fn get(&self, term_id: &HpoTermId) -> Option<&HpoTermInternal> {
        self.hpo_terms.get(term_id)
    }

    pub(crate) fn get_unchecked(&self, term_id: &HpoTermId) -> &HpoTermInternal {
        self.hpo_terms.get_unchecked(term_id)
    }

    fn get_mut(&mut self, term_id: &HpoTermId) -> Option<&mut HpoTermInternal> {
        self.hpo_terms.get_mut(term_id)
    }

    fn get_unchecked_mut(&mut self, term_id: &HpoTermId) -> &mut HpoTermInternal {
        self.hpo_terms.get_unchecked_mut(term_id)
    }

}

/// Methods to add annotations
///
/// These methods should rarely (if ever) be used by clients.
/// Calling these functions might disrupt the Ontology and associated terms.
impl Ontology {
    /// Add a gene to the Ontology. and return the [`GeneId`]
    ///
    /// If the gene does not yet exist, a new [`Gene`] entity is created
    /// and stored in the Ontology.
    /// If the gene already exists in the ontology, it is not added again.
    ///
    /// # Note
    ///
    /// Adding a gene does not connect it to any HPO terms.
    /// Use [`Ontology::link_gene_term`] for creating connections.
    pub fn add_gene(&mut self, gene_name: &str, gene_id: &str) -> OntologyResult<GeneId> {
        let id = GeneId::try_from(gene_id)?;
        match self.genes.entry(id) {
            std::collections::hash_map::Entry::Occupied(_) => Ok(id),
            std::collections::hash_map::Entry::Vacant(entry) => {
                entry.insert(Gene::new(id, gene_name));
                Ok(id)
            }
        }
    }

    /// Add a OMIM disease to the Ontology. and return the [`OmimDiseaseId`]
    ///
    /// If the disease does not yet exist, a new [`OmimDisease`] entity is
    /// created and stored in the Ontology.
    /// If the disease already exists in the ontology, it is not added again.
    /// # Note
    ///
    /// Adding a disease does not connect it to any HPO terms.
    /// Use [`Ontology::link_omim_disease_term`] for creating connections.
    pub fn add_omim_disease(
        &mut self,
        omim_disease_name: &str,
        omim_disease_id: &str,
    ) -> OntologyResult<OmimDiseaseId> {
        let id = OmimDiseaseId::try_from(omim_disease_id)?;
        match self.omim_diseases.entry(id) {
            std::collections::hash_map::Entry::Occupied(_) => Ok(id),
            std::collections::hash_map::Entry::Vacant(entry) => {
                entry.insert(OmimDisease::new(id, omim_disease_name));
                Ok(id)
            }
        }
    }

    /// Add the [`Gene`] as annotation to the [`HpoTerm`]
    ///
    /// The gene will be recursively connected to all parent `HpoTerms` as well.
    ///
    /// This method does not add the HPO-term to the [`Gene`], this must be handled
    /// by the client.
    ///
    /// # Panics
    ///
    /// If the HPO term is not present, the method will panic
    pub fn link_gene_term(&mut self, term_id: &HpoTermId, gene_id: GeneId) {
        let term = self
            .get_mut(term_id)
            .expect("Cannot add gene to non-existing term");

        if term.add_gene(gene_id) {
            // If the gene is already associated to the term, this branch will
            // be skipped. That is desired, because by definition
            // all parent terms are already linked as well
            let parents = term.all_parents().clone();
            for parent in &parents {
                self.link_gene_term(parent, gene_id);
            }
        }
    }

    /// Add the [`OmimDisease`] as annotation to the [`HpoTerm`]
    ///
    /// The disease will be recursively connected to all parent `HpoTerms` as well.
    ///
    /// This method does not add the HPO-term to the [`OmimDisease`], this
    /// must be handled by the client.
    ///
    /// # Panics
    ///
    /// If the HPO term is not present, the method will panic
    pub fn link_omim_disease_term(&mut self, term_id: &HpoTermId, omim_disease_id: OmimDiseaseId) {
        let term = self
            .get_mut(term_id)
            .expect("Cannot add omim_disease to non-existing term");

        if term.add_omim_disease(omim_disease_id) {
            // If the disease is already associated to the term, this branch will
            // be skipped. That is desired, because by definition
            // all parent terms are already linked as well
            let parents = term.all_parents().clone();
            for parent in &parents {
                self.link_omim_disease_term(parent, omim_disease_id);
            }
        }
    }

    /// Calculates the [`crate::InformationContent`]s for every term
    ///
    /// This method should only be called **after** all terms are added,
    /// connected and all genes and diseases are linked as well.
    pub fn calculate_information_content(&mut self) {
        self.calculate_gene_ic();
        self.calculate_omim_disease_ic();
    }

    fn calculate_gene_ic(&mut self) {
        let n_genes = self.genes.len() as f32;

        if n_genes == 0.0 {
            // No genes present in the Ontology
            // so we keep the `Default` value for
            // the Information content
            return;
        }

        for term in self.hpo_terms.values_mut() {
            let ic_gene = match term.genes().len() {
                0 => 0.0,
                n => (n as f32 / n_genes).ln() * -1.0,
            };
            let ic = term.information_content_mut();
            *ic.gene_mut() = ic_gene;
        }
    }

    fn calculate_omim_disease_ic(&mut self) {
        let n_omim_diseases = self.omim_diseases.len() as f32;

        if n_omim_diseases == 0.0 {
            // No omim_diseases present in the Ontology
            // so we keep the `Default` value for
            // the Information content
            return;
        }

        for term in self.hpo_terms.values_mut() {
            let ic_omim_disease = match term.omim_diseases().len() {
                0 => 0.0,
                n => (n as f32 / n_omim_diseases).ln() * -1.0,
            };
            let ic = term.information_content_mut();
            *ic.omim_disease_mut() = ic_omim_disease;
        }
    }
}

/// Crate-only functions for setting up and building the Ontology
///
/// Those methods should not be exposed publicly
impl Ontology {
    fn all_grandparents(&mut self, term_id: &HpoTermId) -> &HpoParents {
        // This looks weird, but I could not find another way to statisfy the Borrow checker
        let cached = {
            let term = self.get_unchecked(term_id);
            term.parents_cached()
        };
        if !cached {
            self.create_cache_of_grandparents(term_id);
        }
        let term = self.get_unchecked(term_id);
        term.all_parents()
    }

    fn create_cache_of_grandparents(&mut self, term_id: &HpoTermId) {
        let term = self.get_unchecked(term_id);
        let parents = term.parents().clone();
        let mut res = HpoParents::default();
        for parent in &parents {
            let grandparents = self.all_grandparents(parent);
            for gp in grandparents {
                res.insert(*gp);
            }
        }
        let term = self.get_unchecked_mut(term_id);
        *term.all_parents_mut() = res.bitor(&parents);
    }

    pub(crate) fn create_cache(&mut self) {
        let term_ids: Vec<HpoTermId> = self.hpo_terms.keys();

        for id in term_ids {
            self.create_cache_of_grandparents(&id);
        }
    }

    pub(crate) fn add_term(&mut self, term: HpoTermInternal) -> HpoTermId {
        let id = *term.id();
        self.hpo_terms.insert(id, term);
        id
    }

    pub(crate) fn add_parent(&mut self, parent_id: HpoTermId, child_id: HpoTermId) {
        let parent = self.get_unchecked_mut(&parent_id);
        parent.add_child(child_id);

        let child = self.get_unchecked_mut(&child_id);
        child.add_parent(parent_id);
    }
}


/// An iterator of [`HpoTerm`]s
pub struct OntologyIterator<'a> {
    inner: std::slice::Iter<'a, HpoTermInternal>,
    ontology: &'a Ontology,
}

impl<'a> std::iter::Iterator for OntologyIterator<'a> {
    type Item = HpoTerm<'a>;
    fn next(&mut self) -> Option<Self::Item> {
        match self.inner.next() {
            Some(term) => Some(HpoTerm::new(self.ontology, term)),
            None => None,
        }
    }
}

impl<'a> IntoIterator for &'a Ontology {
    type Item = HpoTerm<'a>;
    type IntoIter = OntologyIterator<'a>;

    fn into_iter(self) -> Self::IntoIter {
        self.hpos()
    }
}
