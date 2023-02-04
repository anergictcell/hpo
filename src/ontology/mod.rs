use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::Read;
use std::ops::BitOr;
use std::path::Path;

use crate::annotations::{Gene, GeneId};
use crate::annotations::{OmimDisease, OmimDiseaseId};
use crate::parser;
use crate::term::internal::{BinaryTermBuilder, HpoTermInternal};
use crate::term::{HpoParents, HpoTerm};
use crate::u32_from_bytes;
use crate::HpoResult;
use crate::{HpoError, HpoTermId};

use core::fmt::Debug;

mod termarena;
use termarena::Arena;

#[cfg_attr(doc, aquamarine::aquamarine)]
/// `Ontology` is the main interface of the `hpo` crate and contains all data
///
/// The [`Ontology`] struct holds all information about the ontology
/// and the ownership of all [`HpoTerm`]s, [`Gene`]s and [`OmimDisease`]s.
///
/// It is recommended to use the public methods [`Ontology::from_standard`]
/// to build the ontology
/// from standard annotation data from Jax. You will need to download
/// the data from [HPO](https://hpo.jax.org/) itself.
///
/// This crate also provides a snapshot of all relevant data in binary format
/// in `tests/ontology.hpo` which can be loaded via [`Ontology::from_binary`].
/// You should check how up-to-date the snapshot is, though.
///
/// ```mermaid
/// erDiagram
///     ONTOLOGY ||--|{ HPOTERM : contains
///     HPOTERM ||--|{ HPOTERM : is_a
///     HPOTERM }|--o{ DISEASE : phenotype_of
///     HPOTERM }|--o{ GENE : phenotype_of
///     HPOTERM {
///         str name
///         HpoTermId id
///         HpoTerms parents
///         HpoTerms children
///         Genes genes
///         OmimDiseases omim_diseases
///     }
///     DISEASE {
///         str name
///         OmimDiseaseId id
///         HpoGroup hpo_terms
///     }
///     GENE {
///         str name
///         GeneId id
///         HpoGroup hpo_terms
///     }
/// ```
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
    /// - Actual OBO data: [`hp.obo`](https://hpo.jax.org/app/data/ontology)
    /// - Links between HPO and OMIM diseases: [`phenotype.hpoa`](https://hpo.jax.org/app/data/annotations)
    /// - Links between HPO and Genes: [`phenotype_to_genes.txt`](http://purl.obolibrary.org/obo/hp/hpoa/phenotype_to_genes.txt)
    ///
    /// and then specify the folder where the data is stored.
    ///
    /// # Errors
    ///
    /// This method can fail for various reasons:
    ///
    /// - obo file not present or available: [`HpoError::CannotOpenFile`]
    /// - [`Ontology::add_gene`] failed (TODO)
    /// - [`Ontology::add_omim_disease`] failed (TODO)
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use hpo::Ontology;
    /// use hpo::HpoTermId;
    ///
    /// let ontology = Ontology::from_standard("./example_data/").unwrap();
    ///
    /// assert!(ontology.len() > 15_000);
    ///
    /// let absent_term = HpoTermId::try_from("HP:9999999").unwrap();
    /// assert!(ontology.hpo(absent_term).is_none());
    ///
    /// let present_term = HpoTermId::try_from("HP:0000001").unwrap();
    /// let root_term = ontology.hpo(present_term).unwrap();
    /// assert_eq!(root_term.name(), "Phenotypical abnormality");
    /// ```
    ///
    pub fn from_standard(folder: &str) -> HpoResult<Self> {
        let mut ont = Ontology::default();
        let path = Path::new(folder);
        let obo = path.join(crate::OBO_FILENAME);
        let gene = path.join(crate::GENE_FILENAME);
        let disease = path.join(crate::DISEASE_FILENAME);
        parser::load_from_standard_files(&obo, &gene, &disease, &mut ont)?;
        ont.calculate_information_content()?;
        Ok(ont)
    }

    /// Build an Ontology from a binary data blob
    ///
    /// The data must be in the proper format, as defined in
    /// [`Ontology::as_bytes`]. This method adds all terms, creates the
    /// parent-child structure of the ontology, adds genes and Omim diseases
    /// and ensures proper inheritance of gene/disease annotations.
    /// It also calculates the `InformationContent` for every term.
    ///
    /// # Errors
    ///
    /// This method can fail for various reasons:
    ///
    /// - Binary file not available: [`HpoError::CannotOpenFile`]
    /// - `Ontology::add_genes_from_bytes` failed (TODO)
    /// - `Ontology::add_omim_disease_from_bytes` failed (TODO)
    /// - `add_terms_from_bytes` failed (TODO)
    /// - `add_parent_from_bytes` failed (TODO)
    /// - Size of binary data does not match the content: [`HpoError::ParseBinaryError`]
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
    /// assert!(ontology.hpo(absent_term).is_none());
    ///
    /// let present_term = HpoTermId::try_from("HP:0000001").unwrap();
    /// let root_term = ontology.hpo(present_term).unwrap();
    /// assert_eq!(root_term.name(), "All");
    /// ```
    pub fn from_binary<P: AsRef<Path>>(filename: P) -> HpoResult<Self> {
        let mut ont = Ontology::default();
        let bytes = match File::open(filename) {
            Ok(mut file) => {
                let len = file
                    .metadata()
                    .map_err(|_| {
                        HpoError::CannotOpenFile(
                            "unable to get filesize of binary file".to_string(),
                        )
                    })?
                    .len();
                let mut bytes = Vec::with_capacity(len.try_into()?);
                file.read_to_end(&mut bytes).map_err(|_| {
                    HpoError::CannotOpenFile("unable to read from binary file".to_string())
                })?;
                bytes
            }
            Err(_) => {
                return Err(crate::HpoError::CannotOpenFile(
                    "unable to open binary file".to_string(),
                ))
            }
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
            ont.calculate_information_content()?;
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
    ///
    /// # Panics
    ///
    /// Panics when the buffer length of any subsegment larger than `u32::MAX`
    pub fn as_bytes(&self) -> Vec<u8> {
        fn usize_to_u32(n: usize) -> u32 {
            n.try_into().expect("unable to convert {n} to u32")
        }
        let mut res = Vec::new();

        // All HPO Terms
        let mut buffer = Vec::new();
        for term in self.hpo_terms.values() {
            buffer.append(&mut term.as_bytes());
        }
        res.append(&mut usize_to_u32(buffer.len()).to_be_bytes().to_vec());
        res.append(&mut buffer);

        // All Term - Parent connections
        buffer.clear();
        for term in self.hpo_terms.values() {
            buffer.append(&mut term.parents_as_byte());
        }
        res.append(&mut usize_to_u32(buffer.len()).to_be_bytes().to_vec());
        res.append(&mut buffer);

        // Genes and Gene-Term connections
        buffer.clear();
        for gene in self.genes.values() {
            buffer.append(&mut gene.as_bytes());
        }
        res.append(&mut usize_to_u32(buffer.len()).to_be_bytes().to_vec());
        res.append(&mut buffer);

        // OMIM Disease and Disease-Term connections
        buffer.clear();
        for omim_disease in self.omim_diseases.values() {
            buffer.append(&mut omim_disease.as_bytes());
        }
        res.append(&mut usize_to_u32(buffer.len()).to_be_bytes().to_vec());
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
    pub fn hpo(&self, term_id: HpoTermId) -> Option<HpoTerm> {
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

    /// Returns a reference to the [`Gene`] with the provided symbol / name
    ///
    /// If no such gene is present, `None` is returned
    ///
    /// # Note
    ///
    /// `Gene`s are not index by name, so this method searches through all
    /// genes. If you can, prefer using the [`GeneId`] and [`Ontology.gene`].
    pub fn gene_by_name(&self, symbol: &str) -> Option<&Gene> {
        self.genes.values().find(|&gene| gene.name() == symbol)
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

    /// Constructs a smaller ontology that contains only the `leaves` terms and
    /// all terms needed to connect to each leaf to `root`
    ///
    /// # Errors
    ///
    /// Fails if `root` is not an ancestor of all leaves
    pub fn sub_ontology<'a, T:IntoIterator<Item=HpoTerm<'a>>>(&self, root: HpoTerm, leaves: T) -> Result<Self, HpoError> {
        let mut terms = HashSet::new();
        for term in leaves {
            terms.insert(self.get_unchecked(term.id()));
            for parent in term.path_to_ancestor(&root).ok_or(HpoError::NotImplemented)? {
                terms.insert(self.get_unchecked(parent));
            }
        }
        let ids: HashSet<HpoTermId> = terms.iter().map(|term| *term.id()).collect();

        let mut ont = Self::default();
        for term in &terms {
            let internal = HpoTermInternal::new(term.name().to_string(), *term.id());
            ont.add_term(internal);
        };
        for term in &terms {
            for parent in term.parents() {
                if ids.contains(&parent) {
                    ont.add_parent(parent, *term.id());
                }
            }
        };

        ont.create_cache();

        for term in &terms {
            for gene in term.genes() {
                let gene_id = ont.add_gene(
                    self.gene(gene).ok_or(HpoError::DoesNotExist)?.name(),
                    &gene.as_u32().to_string()
                )?;
                ont.link_gene_term(*term.id(), gene_id)?;
                ont
                    .gene_mut(&gene_id).ok_or(HpoError::DoesNotExist)?
                    .add_term(*term.id());
            }

            for omim_disease in term.omim_diseases() {
                let omim_disease_id = ont.add_omim_disease(
                    self.omim_disease(omim_disease).ok_or(HpoError::DoesNotExist)?.name(),
                    &omim_disease.as_u32().to_string()
                )?;
                ont.link_omim_disease_term(*term.id(), omim_disease_id)?;
                ont
                    .omim_disease_mut(&omim_disease_id).ok_or(HpoError::DoesNotExist)?
                    .add_term(*term.id());
            }
        }
        ont.calculate_information_content()?;

        Ok(ont)
    }

    /// Returns the code to crate a `Mermaid` flow diagram
    ///
    /// This is meant to be used with smaller ontologies, e.g. from [`Ontology::sub_ontology`]
    pub fn as_mermaid(&self) -> String {
        let mut code = String::new();
        code.push_str("graph TD\n");
        for term in self {
            code.push_str(&format!("{}[\"{}\n{}\"]\n", term.id(), term.id(), term.name()));
            for child in term.children(){
                code.push_str(&format!("{} --> {}\n", term.id(), child.id()));
            }
        }
        code
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
    ///
    /// # Errors
    ///
    /// If the `gene_id` is invalid, an [`HpoError::ParseIntError`] is returned
    pub fn add_gene(&mut self, gene_name: &str, gene_id: &str) -> HpoResult<GeneId> {
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
    ///
    /// # Note
    ///
    /// Adding a disease does not connect it to any HPO terms.
    /// Use [`Ontology::link_omim_disease_term`] for creating connections.
    ///
    /// # Errors
    ///
    /// If the `omim_disease_id` is invalid, an [`HpoError::ParseIntError`] is returned
    pub fn add_omim_disease(
        &mut self,
        omim_disease_name: &str,
        omim_disease_id: &str,
    ) -> HpoResult<OmimDiseaseId> {
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
    /// # Errors
    ///
    /// If the HPO term is not present, an [`HpoError::DoesNotExist`] is returned
    pub fn link_gene_term(&mut self, term_id: HpoTermId, gene_id: GeneId) -> HpoResult<()> {
        let term = self.get_mut(term_id).ok_or(HpoError::DoesNotExist)?;

        if term.add_gene(gene_id) {
            // If the gene is already associated to the term, this branch will
            // be skipped. That is desired, because by definition
            // all parent terms are already linked as well
            let parents = term.all_parents().clone();
            for parent in &parents {
                self.link_gene_term(parent, gene_id)?;
            }
        }
        Ok(())
    }

    /// Add the [`OmimDisease`] as annotation to the [`HpoTerm`]
    ///
    /// The disease will be recursively connected to all parent `HpoTerms` as well.
    ///
    /// This method does not add the HPO-term to the [`OmimDisease`], this
    /// must be handled by the client.
    ///
    /// # Errors
    ///
    /// If the HPO term is not present, an [`HpoError`] is returned
    pub fn link_omim_disease_term(
        &mut self,
        term_id: HpoTermId,
        omim_disease_id: OmimDiseaseId,
    ) -> HpoResult<()> {
        let term = self.get_mut(term_id).ok_or(HpoError::DoesNotExist)?;

        if term.add_omim_disease(omim_disease_id) {
            // If the disease is already associated to the term, this branch will
            // be skipped. That is desired, because by definition
            // all parent terms are already linked as well
            let parents = term.all_parents().clone();
            for parent in &parents {
                self.link_omim_disease_term(parent, omim_disease_id)?;
            }
        }
        Ok(())
    }

    /// Calculates the [`crate::term::InformationContent`]s for every term
    ///
    /// This method should only be called **after** all terms are added,
    /// connected and all genes and diseases are linked as well.
    ///
    /// It can be called repeatedly, all values are recalculated each time,
    /// as long as the Ontology contains at least 1 gene/disease.
    /// When no genes/diseases are present, the IC is not calculated nor updated.
    ///
    /// # Errors
    ///
    /// This method returns an error if there are more Genes or Terms than `u16::MAX`
    /// because larger numbers can't be safely converted to `f32`
    pub fn calculate_information_content(&mut self) -> HpoResult<()> {
        self.calculate_gene_ic()?;
        self.calculate_omim_disease_ic()?;
        Ok(())
    }

    /// Calculates the gene-specific Information Content for every term
    ///
    /// If no genes are present in the Ontology, no IC are calculated
    fn calculate_gene_ic(&mut self) -> HpoResult<()> {
        let n_genes = self.genes.len();
        for term in self.hpo_terms.values_mut() {
            let current_genes = term.genes().len();
            term.information_content_mut()
                .set_gene(n_genes, current_genes)?;
        }
        Ok(())
    }

    /// Calculates the Omim-Disease-specific Information Content for every term
    ///
    /// If no diseases are present in the Ontology, no IC are calculated
    fn calculate_omim_disease_ic(&mut self) -> HpoResult<()> {
        let n_omim_diseases = self.omim_diseases.len();

        for term in self.hpo_terms.values_mut() {
            let current_diseases = term.omim_diseases().len();
            term.information_content_mut()
                .set_omim_disease(n_omim_diseases, current_diseases)?;
        }
        Ok(())
    }
}

/// Crate-only functions for setting up and building the Ontology
///
/// Those methods should not be exposed publicly
impl Ontology {
    /// Adds an [`HpoTerm`] to the ontology
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

    /// Connects an [`HpoTerm`] to its parent term
    ///
    /// This method is part of the Ontology-building, based on the binary
    /// data format and requires a specified data layout.
    ///
    /// The method assumes that the data is in the right format and also
    /// assumes that the caller will populate the `all_parents` caches for
    /// each term.
    ///
    /// See [`HpoTermInternal::parents_as_byte`] for explanation of the binary layout.
    ///
    /// # Panics
    ///
    /// This method will panic if the length of bytes does not exactly correspond
    /// to the contained data
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
    fn add_genes_from_bytes(&mut self, bytes: &[u8]) -> HpoResult<()> {
        let mut idx: usize = 0;
        loop {
            if idx >= bytes.len() {
                break;
            }
            let gene_len = u32_from_bytes(&bytes[idx..]) as usize;
            let gene = Gene::try_from(&bytes[idx..idx + gene_len])?;
            for term in gene.hpo_terms() {
                self.link_gene_term(term, *gene.id())?;
            }
            self.genes.insert(*gene.id(), gene);
            idx += gene_len;
        }
        Ok(())
    }

    /// Adds [`OmimDisease`]s to the ontoloigy and connects them to connected terms
    ///
    /// This method is part of the Ontology-building, based on the binary
    /// data format and requires a specified data layout.
    ///
    /// It connects all connected terms and their parents properly. The
    /// method assumes that the bytes encode all Disease-term connections.
    ///
    /// See [`OmimDisease::as_bytes`] for explanation of the binary layout
    fn add_omim_disease_from_bytes(&mut self, bytes: &[u8]) -> HpoResult<()> {
        let mut idx: usize = 0;
        loop {
            if idx >= bytes.len() {
                break;
            }
            let disease_len = u32_from_bytes(&bytes[idx..]) as usize;
            let disease = OmimDisease::try_from(&bytes[idx..idx + disease_len])?;
            for term in disease.hpo_terms() {
                self.link_omim_disease_term(term, *disease.id())?;
            }
            self.omim_diseases.insert(*disease.id(), disease);
            idx += disease_len;
        }
        Ok(())
    }

    /// This method is part of the cache creation to link all terms to their
    /// direct and indirect parents (grandparents)
    ///
    /// # Panics
    ///
    /// This method will panic if the `term_id` is not present in the Ontology
    fn all_grandparents(&mut self, term_id: HpoTermId) -> &HpoParents {
        if !self.get_unchecked(term_id).parents_cached() {
            self.create_cache_of_grandparents(term_id);
        }
        let term = self.get_unchecked(term_id);
        term.all_parents()
    }

    /// This method is part of the cache creation to link all terms to their
    /// direct and indirect parents (grandparents)
    ///
    /// It will (somewhat) recursively iterate all parents and copy all their parents.
    /// During this recursion, the list of `all_parents` is cached in each term that was
    /// iterated.
    ///
    /// The logic is that the recursion bubbles up all the way to the top of the ontolgy
    /// and then caches the list of direct and indirect parents for every term bubbling
    /// back down. The recursion does not reach the top level again, because it will stop
    /// once it reaches a term with already cached `all_parents`.
    ///
    /// # Panics
    ///
    /// This method will panic if the `term_id` is not present in the Ontology
    fn create_cache_of_grandparents(&mut self, term_id: HpoTermId) {
        let mut res = HpoParents::default();
        let parents = self.get_unchecked(term_id).parents().clone();
        for parent in &parents {
            let grandparents = self.all_grandparents(parent);
            for gp in grandparents {
                res.insert(gp);
            }
        }
        let term = self.get_unchecked_mut(term_id);
        *term.all_parents_mut() = res.bitor(&parents);
    }

    /// Crates and caches the `all_parents` values for every term
    ///
    /// This method can only be called once and afterwards no new terms
    /// should be added to the Ontology anymore and no new term-parent connection
    /// should be created.
    /// Since this method caches the results, rerunning it will not cause a new
    /// calculation.
    pub(crate) fn create_cache(&mut self) {
        let term_ids: Vec<HpoTermId> = self.hpo_terms.keys();

        for id in term_ids {
            self.create_cache_of_grandparents(id);
        }
    }

    /// Insert an `HpoTermInternal` to the ontology
    ///
    /// This method does not link the term to its parents or to any annotations
    pub(crate) fn add_term(&mut self, term: HpoTermInternal) -> HpoTermId {
        let id = *term.id();
        self.hpo_terms.insert(term);
        id
    }

    /// Add a connection from an [`HpoTerm`] to its parent
    ///
    /// This method is called once for every dependency in the Ontology during the initialization.
    ///
    /// There should rarely be a need to call this method outside of the ontology building
    ///
    /// # Panics
    ///
    /// This method will panic if the `parent_id` or `child_id` is not present in the Ontology
    pub(crate) fn add_parent(&mut self, parent_id: HpoTermId, child_id: HpoTermId) {
        let parent = self.get_unchecked_mut(parent_id);
        parent.add_child(child_id);

        let child = self.get_unchecked_mut(child_id);
        child.add_parent(parent_id);
    }

    /// Returns the `HpoTermInternal` with the given `HpoTermId`
    ///
    /// Returns `None` if no such term is present
    pub(crate) fn get(&self, term_id: HpoTermId) -> Option<&HpoTermInternal> {
        self.hpo_terms.get(term_id)
    }

    /// Returns the `HpoTermInternal` with the given `HpoTermId`
    ///
    /// This method should only be called if the caller is sure that the term actually
    /// exists, e.g. during an iteration of all `HpoTermId`s.
    ///
    /// # Panics
    ///
    /// This method will panic if the `term_id` is not present in the Ontology
    pub(crate) fn get_unchecked(&self, term_id: HpoTermId) -> &HpoTermInternal {
        self.hpo_terms.get_unchecked(term_id)
    }

    /// Returns a mutable reference to the `HpoTermInternal` with the given `HpoTermId`
    ///
    /// Returns `None` if no such term is present
    fn get_mut(&mut self, term_id: HpoTermId) -> Option<&mut HpoTermInternal> {
        self.hpo_terms.get_mut(term_id)
    }

    /// Returns a mutable reference to the `HpoTermInternal` with the given `HpoTermId`
    ///
    /// This method should only be called if the caller is sure that the term actually
    /// exists, e.g. during an iteration of all `HpoTermId`s.
    ///
    /// # Panics
    ///
    /// This method will panic if the `term_id` is not present in the Ontology
    fn get_unchecked_mut(&mut self, term_id: HpoTermId) -> &mut HpoTermInternal {
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

        assert_eq!(ont.get_unchecked(3u32.into()).parents().len(), 2);
        assert_eq!(ont.get_unchecked(1u32.into()).children().len(), 1);
        assert_eq!(ont.get_unchecked(2u32.into()).children().len(), 1);
    }
}
