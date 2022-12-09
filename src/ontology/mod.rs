use std::collections::HashMap;
use std::ops::BitOr;

use crate::annotations::{Gene, GeneId};
use crate::annotations::{OmimDisease, OmimDiseaseId};
use crate::term::HpoTerm;
use crate::term::HpoTermInternal;
use crate::HpoParents;
use crate::HpoTermId;
use crate::OntologyResult;

use core::fmt::Debug;

mod termarena;
use termarena::Arena;

#[derive(Default)]
pub struct Ontology {
    hpo_terms: Arena,
    hpo_ids: Vec<HpoTermId>,
    genes: HashMap<GeneId, Gene>,
    omim_diseases: HashMap<OmimDiseaseId, OmimDisease>,
}

impl Debug for Ontology {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Ontology with {} terns", self.hpo_terms.len())
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

    pub fn create_cache(&mut self) {
        let term_ids: Vec<HpoTermId> = self.hpo_terms.keys();

        for id in term_ids {
            self.create_cache_of_grandparents(&id);
        }
    }

    fn add_term(&mut self, term: HpoTermInternal) -> HpoTermId {
        let id = *term.id();
        self.hpo_terms.insert(id, term);
        self.hpo_ids.push(id);
        id
    }

    pub fn add_term_by_name(&mut self, name: &str) -> HpoTermId {
        let t = HpoTermInternal::new(name);
        self.add_term(t)
    }

    pub fn add_parent(&mut self, parent_id: HpoTermId, child_id: HpoTermId) {
        let parent = self.get_unchecked_mut(&parent_id);
        parent.add_child(child_id);

        let child = self.get_unchecked_mut(&child_id);
        child.add_parent(parent_id);
    }

    /// Is this method really needed or beneficial?
    pub fn shrink_to_fit(&mut self) {
        self.hpo_terms.shrink_to_fit();
    }

    pub fn len(&self) -> usize {
        self.hpo_terms.len()
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

/// Methods to add annotations
impl Ontology {
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

    pub fn link_gene_term(&mut self, term_id: &HpoTermId, gene_id: GeneId) {
        let term = self
            .get_mut(term_id)
            .expect("Cannot add gene to non-existing term");

        if term.add_gene(gene_id) {
            // this part can be skipped, if the gene is already linked to the term,
            // because all parent terms are already linked as well
            let parents = term.all_parents().clone();
            for parent in &parents {
                self.link_gene_term(parent, gene_id);
            }
        }
    }

    pub fn link_omim_disease_term(&mut self, term_id: &HpoTermId, omim_disease_id: OmimDiseaseId) {
        let term = self
            .get_mut(term_id)
            .expect("Cannot add omim_disease to non-existing term");

        if term.add_omim_disease(omim_disease_id) {
            // this part can be skipped, if the omim_disease is already linked to the term,
            // because all parent terms are already linked as well
            let parents = term.all_parents().clone();
            for parent in &parents {
                self.link_omim_disease_term(parent, omim_disease_id);
            }
        }
    }

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

/// Public API of the Ontology
///
/// Those methods are all safe to use
impl Ontology {
    pub fn get(&self, term_id: &HpoTermId) -> Option<&HpoTermInternal> {
        self.hpo_terms.get(term_id)
    }

    pub fn get_unchecked(&self, term_id: &HpoTermId) -> &HpoTermInternal {
        self.hpo_terms.get_unchecked(term_id)
    }

    fn get_mut(&mut self, term_id: &HpoTermId) -> Option<&mut HpoTermInternal> {
        self.hpo_terms.get_mut(term_id)
    }

    fn get_unchecked_mut(&mut self, term_id: &HpoTermId) -> &mut HpoTermInternal {
        self.hpo_terms.get_unchecked_mut(term_id)
    }

    pub fn hpo(&self, term_id: &HpoTermId) -> OntologyResult<HpoTerm> {
        HpoTerm::try_new(self, term_id)
    }

    pub fn hpos(&self) -> OntologyIterator {
        OntologyIterator {
            inner: self.hpo_terms.values().iter(),
            ontology: self,
        }
    }

    pub fn gene(&self, gene_id: &GeneId) -> Option<&Gene> {
        self.genes.get(gene_id)
    }

    pub fn gene_mut(&mut self, gene_id: &GeneId) -> Option<&mut Gene> {
        self.genes.get_mut(gene_id)
    }

    pub fn genes(&self) -> std::collections::hash_map::Values<'_, GeneId, Gene> {
        self.genes.values()
    }

    pub fn omim_disease(&self, omim_disease_id: &OmimDiseaseId) -> Option<&OmimDisease> {
        self.omim_diseases.get(omim_disease_id)
    }

    pub fn omim_disease_mut(
        &mut self,
        omim_disease_id: &OmimDiseaseId,
    ) -> Option<&mut OmimDisease> {
        self.omim_diseases.get_mut(omim_disease_id)
    }

    pub fn omim_diseases(
        &self,
    ) -> std::collections::hash_map::Values<'_, OmimDiseaseId, OmimDisease> {
        self.omim_diseases.values()
    }
}

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
