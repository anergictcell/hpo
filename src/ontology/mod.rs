//! Ontology is the heart and main interface of the full crate.
//!
//! The [`Ontology`] struct holds all information about the ontology
//! and the ownership of all [`HpoTerm`]s, [`Gene`]s and [`OmimDisease`]s.

use std::collections::HashMap;
use std::fs::File;
use std::io::Read;
use std::ops::BitOr;
use std::path::Path;

use crate::annotations::{Gene, GeneId};
use crate::annotations::{OmimDisease, OmimDiseaseId};
use crate::term::internal::{BinaryTermBuilder, HpoTermInternal};
use crate::term::HpoTerm;
use crate::u32_from_bytes;
use crate::OntologyResult;
use crate::{parser, HpoParents};
use crate::{HpoError, HpoTermId};

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
        ont.calculate_information_content();
        ont
    }

    /// Build an Ontology from a binary data blob
    ///
    /// The data must be in the proper format, as defined in
    /// [`Ontology::as_bytes`]. This method adds all terms, creates the
    /// parent-child structure of the ontology, adds genes and Omim diseases
    /// and ensures proper inheritance of gene/disease annotations.
    /// It also calculates the InformationContent for every term.
    ///
    /// # Examples
    ///
    /// ```
    /// use hpo::{Ontology, HpoTermId};
    ///
    /// let ontology = Ontology::from_binary("./tests/ontology.hpo").unwrap();
    ///
    /// assert!(ontology.len() > 15_000);
    ///
    /// let absent_term = HpoTermId::try_from("HP:9999999").unwrap();
    /// assert!(ontology.hpo(&absent_term).is_none());
    ///
    /// let present_term = HpoTermId::try_from("HP:0000001").unwrap();
    /// let root_term = ontology.hpo(&present_term).unwrap();
    /// assert_eq!(root_term.name(), "All");
    /// ```
    pub fn from_binary<P: AsRef<Path>>(filename: P) -> OntologyResult<Self> {
        let mut ont = Ontology::default();
        let bytes = match File::open(filename) {
            Ok(mut file) => {
                let len = file.metadata().map_err(|_| HpoError::DoesNotExist)?.len();
                let mut bytes = Vec::with_capacity(len.try_into().unwrap());
                file.read_to_end(&mut bytes)
                    .map_err(|_| HpoError::DoesNotExist)?;
                bytes
            }
            Err(_) => return Err(crate::HpoError::DoesNotExist),
        };

        let mut section_start = 0;
        let mut section_end: usize;

        // Terms
        let mut section_len = u32_from_bytes(&bytes[section_start..]) as usize;
        section_end = 4 + section_len;
        ont.add_terms_from_bytes(&bytes[4..section_end]);
        section_start += section_len + 4;

        // Term - Parents
        section_len = u32_from_bytes(&bytes[section_start..]) as usize;
        section_end += 4 + section_len;
        ont.add_parent_from_bytes(&bytes[section_start + 4..section_end]);
        ont.create_cache();
        section_start += section_len + 4;

        // Genes
        section_len = u32_from_bytes(&bytes[section_start..]) as usize;
        section_end += 4 + section_len;
        ont.add_genes_from_bytes(&bytes[section_start + 4..section_end])?;
        section_start += section_len + 4;

        // Omim Diseases
        section_len = u32_from_bytes(&bytes[section_start..]) as usize;
        section_end += 4 + section_len;
        ont.add_omim_disease_from_bytes(&bytes[section_start + 4..section_end])?;
        section_start += section_len + 4;

        if section_start == bytes.len() {
            ont.calculate_information_content();
            Ok(ont)
        } else {
            Err(HpoError::ParseBinaryError)
        }
    }

    /// Returns a binary representation of the Ontology
    ///
    /// The binary data is separated into sections:
    ///
    /// - Terms (Names + IDs) (see `HpoTermInternal::as_bytes`)
    /// - Term - Parent connection (Child ID - Parent ID)
    ///   (see `HpoTermInternal::parents_as_byte`)
    /// - Genes (Names + IDs + Connected HPO Terms) ([`Gene::as_bytes`])
    /// - OMIM Diseases (Names + IDs + Connected HPO Terms)
    ///   ([`OmimDisease::as_bytes`])
    ///
    /// Every section starts with 4 bytes to indicate its size
    /// (big-endian encoded `u32`)
    ///
    /// This method is only useful if you use are modifying the ontology
    /// and want to save data for later re-use.
    pub fn as_bytes(&self) -> Vec<u8> {
        let mut res = Vec::new();

        // All HPO Terms
        let mut buffer = Vec::new();
        for term in self.hpo_terms.values() {
            buffer.append(&mut term.as_bytes());
        }
        res.append(&mut (buffer.len() as u32).to_be_bytes().to_vec());
        res.append(&mut buffer);

        // All Term - Parent connections
        buffer.clear();
        for term in self.hpo_terms.values() {
            buffer.append(&mut term.parents_as_byte());
        }
        res.append(&mut (buffer.len() as u32).to_be_bytes().to_vec());
        res.append(&mut buffer);

        // Genes and Gene-Term connections
        buffer.clear();
        for gene in self.genes.values() {
            buffer.append(&mut gene.as_bytes());
        }
        res.append(&mut (buffer.len() as u32).to_be_bytes().to_vec());
        res.append(&mut buffer);

        // OMIM Disease and Disease-Term connections
        buffer.clear();
        for omim_disease in self.omim_diseases.values() {
            buffer.append(&mut omim_disease.as_bytes());
        }
        res.append(&mut (buffer.len() as u32).to_be_bytes().to_vec());
        res.append(&mut buffer);

        res
    }

    /// Returns the number of HPO-Terms in the Ontology
    pub fn len(&self) -> usize {
        self.hpo_terms.len()
    }

    /// Returns `true` if the Ontology does not contain any HPO-Terms
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

    /// Returns a reference to the [`Gene`] of the provided [`GeneId`]
    ///
    /// If no such gene is present, `None` is returned
    pub fn gene(&self, gene_id: &GeneId) -> Option<&Gene> {
        self.genes.get(gene_id)
    }

    /// Returns a mutable reference to the [`Gene`] of the provided [`GeneId`]
    ///
    /// If no such gene is present, `None` is returned
    pub fn gene_mut(&mut self, gene_id: &GeneId) -> Option<&mut Gene> {
        self.genes.get_mut(gene_id)
    }

    /// Returns an Iterator of all [`Gene`]s from the Ontology
    ///
    /// It is likely that the return type will change to a dedicated Iterator
    pub fn genes(&self) -> std::collections::hash_map::Values<'_, GeneId, Gene> {
        self.genes.values()
    }

    /// Returns a reference to the [`OmimDisease`] of the provided [`OmimDiseaseId`]
    ///
    /// If no such disease is present, `None` is returned
    pub fn omim_disease(&self, omim_disease_id: &OmimDiseaseId) -> Option<&OmimDisease> {
        self.omim_diseases.get(omim_disease_id)
    }

    /// Returns a mutable reference to the [`OmimDisease`] of the provided [`OmimDiseaseId`]
    ///
    /// If no such disease is present, `None` is returned
    pub fn omim_disease_mut(
        &mut self,
        omim_disease_id: &OmimDiseaseId,
    ) -> Option<&mut OmimDisease> {
        self.omim_diseases.get_mut(omim_disease_id)
    }

    /// Returns an Iterator of all [`OmimDisease`]s from the Ontology
    ///
    /// It is likely that the return type will change to a dedicated Iterator
    pub fn omim_diseases(
        &self,
    ) -> std::collections::hash_map::Values<'_, OmimDiseaseId, OmimDisease> {
        self.omim_diseases.values()
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
    /// Adds an HpoTerm to the ontology
    ///
    /// This method is part of the Ontology-building, based on the binary
    /// data format and requires a specified data layout.
    ///
    /// The method assumes that the data is in the right format and also
    /// assumes that the caller takes care of handling all consistencies
    /// like parent-child connection etc.
    ///
    /// See [`HpoTermInternal::as_bytes`] for explanation of the binary layout.
    fn add_terms_from_bytes(&mut self, bytes: &[u8]) {
        for term in BinaryTermBuilder::new(bytes) {
            self.add_term(term);
        }
    }

    /// Connects an HpoTerm to its parent term
    ///
    /// This method is part of the Ontology-building, based on the binary
    /// data format and requires a specified data layout.
    ///
    /// The method assumes that the data is in the right format and also
    /// assumes that the caller will populate the all_parents caches for
    /// each term.
    ///
    /// See [`HpoTermInternal::parents_as_byte`] for explanation of the binary layout.
    fn add_parent_from_bytes(&mut self, bytes: &[u8]) {
        let mut idx: usize = 0;
        loop {
            if idx == bytes.len() {
                break;
            }
            let n_parents = u32_from_bytes(&bytes[idx..]) as usize;
            idx += 4;
            let term =
                HpoTermId::from([bytes[idx], bytes[idx + 1], bytes[idx + 2], bytes[idx + 3]]);
            idx += 4;
            for _ in 0..n_parents {
                let parent =
                    HpoTermId::from([bytes[idx], bytes[idx + 1], bytes[idx + 2], bytes[idx + 3]]);
                self.add_parent(parent, term);
                idx += 4;
            }
        }
    }

    /// Adds genes to the ontoloigy and connects them to connected terms
    ///
    /// This method is part of the Ontology-building, based on the binary
    /// data format and requires a specified data layout.
    ///
    /// It connects all connected terms and their parents properly. The
    /// method assumes that the bytes encode all gene-term connections.
    ///
    /// See [`Gene::as_bytes`] for explanation of the binary layout
    fn add_genes_from_bytes(&mut self, bytes: &[u8]) -> OntologyResult<()> {
        let mut idx: usize = 0;
        loop {
            if idx >= bytes.len() {
                break;
            }
            let gene_len = u32_from_bytes(&bytes[idx..]) as usize;
            let gene = Gene::try_from(&bytes[idx..idx + gene_len])?;
            for term in gene.hpo_terms() {
                self.link_gene_term(term, *gene.id())
            }
            self.genes.insert(*gene.id(), gene);
            idx += gene_len;
        }
        Ok(())
    }

    /// Adds OmimDiseases to the ontoloigy and connects them to connected terms
    ///
    /// This method is part of the Ontology-building, based on the binary
    /// data format and requires a specified data layout.
    ///
    /// It connects all connected terms and their parents properly. The
    /// method assumes that the bytes encode all Disease-term connections.
    ///
    /// See [`OmimDisease::as_bytes`] for explanation of the binary layout
    fn add_omim_disease_from_bytes(&mut self, bytes: &[u8]) -> OntologyResult<()> {
        let mut idx: usize = 0;
        loop {
            if idx >= bytes.len() {
                break;
            }
            let disease_len = u32_from_bytes(&bytes[idx..]) as usize;
            let disease = OmimDisease::try_from(&bytes[idx..idx + disease_len])?;
            for term in disease.hpo_terms() {
                self.link_omim_disease_term(term, *disease.id())
            }
            self.omim_diseases.insert(*disease.id(), disease);
            idx += disease_len;
        }
        Ok(())
    }

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

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn add_terms() {
        let test_terms = [
            ("t1", 1u32),
            ("Term with a very long name", 2u32),
            ("", 3u32),
            ("Abnormality", 4u32),
        ];

        let mut ont = Ontology::default();

        let mut v: Vec<u8> = Vec::new();
        for (name, id) in test_terms {
            let t = HpoTermInternal::new(String::from(name), id.into());
            v.append(&mut t.as_bytes());
        }
        ont.add_terms_from_bytes(&v);
        assert_eq!(ont.len(), 4);
    }

    #[test]
    fn add_parents() {
        let test_terms = [
            ("t1", 1u32),
            ("Term with a very long name", 2u32),
            ("", 3u32),
            ("Abnormality", 4u32),
        ];

        let mut ont = Ontology::default();

        let mut v: Vec<u8> = Vec::new();
        for (name, id) in test_terms {
            let t = HpoTermInternal::new(String::from(name), id.into());
            v.append(&mut t.as_bytes());
        }
        ont.add_terms_from_bytes(&v);
        assert_eq!(ont.len(), 4);

        // The fake term has the same HpoTermId as one of of the Test ontology
        let mut fake_term = HpoTermInternal::new(String::from(""), 3u32.into());
        fake_term.add_parent(1u32.into());
        fake_term.add_parent(2u32.into());

        let bytes = fake_term.parents_as_byte();

        ont.add_parent_from_bytes(&bytes[..]);

        assert_eq!(ont.get_unchecked(&3u32.into()).parents().len(), 2);
        assert_eq!(ont.get_unchecked(&1u32.into()).children().len(), 1);
        assert_eq!(ont.get_unchecked(&2u32.into()).children().len(), 1);
    }
}
