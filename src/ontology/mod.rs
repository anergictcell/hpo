use std::collections::{HashMap, HashSet};
use std::ops::BitOr;
use std::slice::Iter;

use crate::term::HpoParentIterator;
use crate::term::HpoTerm;
use crate::term::HpoTermInternal;
use crate::HpoParents;
use crate::HpoTermId;
use crate::OntologyResult;
use crate::annotations::{GeneId, Gene};

use core::fmt::Debug;

#[derive(Default)]
pub struct Ontology {
    hpo_terms: HashMap<HpoTermId, HpoTermInternal>,
    hpo_ids: Vec<HpoTermId>,
    genes: HashMap<GeneId, Gene>
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
            let term = self.get(term_id).unwrap();
            term.parents_cached()
        };
        if !cached {
            self.create_cache_of_grandparents(term_id);
        }
        let term = self.get(term_id).unwrap();
        term.all_parents()
    }

    fn create_cache_of_grandparents(&mut self, term_id: &HpoTermId) {
        let term = self.get(term_id).unwrap();
        let parents = term.parents().clone();
        let mut res = HashSet::default();
        for parent in &parents {
            let grandparents = self.all_grandparents(parent);
            for gp in grandparents {
                res.insert(*gp);
            }
        }
        let term = self.get_mut(term_id).unwrap();
        *term.all_parents_mut() = res.bitor(&parents);
    }

    pub fn create_cache(&mut self) {
        let term_ids: Vec<HpoTermId> = self.hpo_terms.keys().copied().collect();

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
        match self.hpo_terms.get_mut(&parent_id) {
            Some(term) => term.add_child(child_id),
            None => panic!("term not present"),
        }

        match self.hpo_terms.get_mut(&child_id) {
            Some(term) => term.add_parent(parent_id),
            None => panic!("term not present"),
        }
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

    pub fn link_gene_term(&mut self, term_id: &HpoTermId, gene_id: GeneId) {
        let term = self.get_mut(term_id).expect("Cannot add gene to non-existing term");

        // if the gene is already linked to the term, all parent terms are already linked as well
        if term.add_gene(gene_id) {
            let parents = term.all_parents().clone();
            for parent in parents {
                self.link_gene_term(&parent, gene_id);
            }
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

    fn get_mut(&mut self, term_id: &HpoTermId) -> Option<&mut HpoTermInternal> {
        self.hpo_terms.get_mut(term_id)
    }

    pub fn get_term(&self, term_id: &HpoTermId) -> OntologyResult<HpoTerm> {
        HpoTerm::try_new(self, term_id)
    }

    pub fn get_gene(&self, gene_id: &GeneId) -> Option<&Gene> {
        self.genes.get(gene_id)
    }

    /// Returns an iterator over all direct parents of the term
    ///
    /// Same as `HpoTerm.parents()`, but (one would think) it is more performant.
    /// Turns out, there is no measurable difference
    pub fn parents(&self, term_id: &HpoTermId) -> HpoParentIterator {
        let term = self.get(term_id).unwrap();
        HpoParentIterator::new(term.parents(), self)
    }

    /// Returns an iterator over all direct parents of the term
    ///
    /// Same as `HpoTerm.parents()`, but (one would think) it is more performant.
    /// Turns out, there is no measurable difference
    pub fn all_parents(&self, term_id: &HpoTermId) -> HpoParentIterator {
        let term = self.get(term_id).unwrap();
        HpoParentIterator::new(term.all_parents(), self)
    }

    pub fn common_ancestors(&self, t1: &HpoTermId, t2: &HpoTermId) -> HpoParents {
        let term1 = self.get(t1).unwrap();
        let term2 = self.get(t2).unwrap();
        term1.all_parents().bitor(term2.all_parents())
    }

    pub fn iter(&self) -> Iter<HpoTermId> {
        self.hpo_ids.iter()
    }
}
