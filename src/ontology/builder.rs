use crate::annotations::Disease;
use crate::term::internal::HpoTermInternal;
use std::collections::hash_map::Entry;
use std::collections::HashMap;
use std::marker::PhantomData;
use std::ops::BitOr;

use crate::annotations::{Gene, GeneId};
use crate::annotations::{OmimDisease, OmimDiseaseId};
use crate::annotations::{OrphaDisease, OrphaDiseaseId};
use crate::parser::binary::{BinaryTermBuilder, BinaryVersion, Bytes};
use crate::term::HpoGroup;

use crate::{u32_from_bytes, HpoTermId, Ontology};
use crate::HpoResult;
use crate::HpoError;

use crate::ontology::termarena::Arena;


pub struct LooseCollection;
pub struct AllTerms;
pub struct ConnectedTerms;
pub struct FullyAnnotated;

pub trait AddAnotation{}
impl AddAnotation for LooseCollection{}
impl AddAnotation for AllTerms{}
impl AddAnotation for ConnectedTerms{}

fn transition_state<TX, TY>(builder: Builder<TX>) -> Builder<TY> {
    Builder::<TY>{
        hpo_terms: builder.hpo_terms,
        genes: builder.genes,
        omim_diseases: builder.omim_diseases,
        orpha_diseases: builder.orpha_diseases,
        hpo_version: builder.hpo_version,
        categories: builder.categories,
        modifier: builder.modifier,
        state: PhantomData
    }
}


pub struct Builder<T> {
    hpo_terms: Arena,
    genes: HashMap<GeneId, Gene>,
    omim_diseases: HashMap<OmimDiseaseId, OmimDisease>,
    orpha_diseases: HashMap<OrphaDiseaseId, OrphaDisease>,
    hpo_version: (u16, u8, u8),
    categories: HpoGroup,
    modifier: HpoGroup,
    state: PhantomData<T>
}


impl<T: AddAnotation> Builder<T> {
    pub fn add_gene(&mut self, gene_name: &str, gene_id: GeneId) {
        if let Entry::Vacant(entry) = self.genes.entry(gene_id) {
            entry.insert(Gene::new(gene_id, gene_name));
        }
    }

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

    pub fn add_orpha_disease(
        &mut self,
        orpha_disease_name: &str,
        orpha_disease_id: &str,
    ) -> HpoResult<OrphaDiseaseId> {
        let id = OrphaDiseaseId::try_from(orpha_disease_id)?;
        match self.orpha_diseases.entry(id) {
            std::collections::hash_map::Entry::Occupied(_) => Ok(id),
            std::collections::hash_map::Entry::Vacant(entry) => {
                entry.insert(OrphaDisease::new(id, orpha_disease_name));
                Ok(id)
            }
        }
    }
}



impl Builder<LooseCollection>{
    pub fn new() -> Builder<LooseCollection> {
        Builder::<LooseCollection> {
            hpo_terms: Arena::default(),
            genes: HashMap::default(),
            omim_diseases: HashMap::default(),
            orpha_diseases: HashMap::default(),
            hpo_version: (0u16, 0u8, 0u8),
            categories: HpoGroup::default(),
            modifier: HpoGroup::default(),
            state: PhantomData
        }
    }    

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
    pub (crate) fn add_terms_from_bytes(&mut self, bytes: Bytes) {
        for term in BinaryTermBuilder::new(bytes) {
            self.add_term(term);
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

    #[must_use]
    pub fn terms_complete(self) -> Builder<AllTerms> {
        transition_state(self)
    }
}

impl Builder<AllTerms> {
    /// Add a connection from an [`HpoTerm`] to its parent
    ///
    /// This method is called once for every dependency in the Ontology during the initialization.
    ///
    /// There should rarely be a need to call this method outside of the ontology building
    ///
    /// # Panics
    ///
    /// This method will panic if the `parent_id` or `child_id` is not present in the Ontology
    ///
    /// # Examples
    ///
    /// ```
    /// use hpo::Ontology;
    ///
    /// let mut ontology = Ontology::default();
    /// ontology.insert_term("Foo".into(), 1u32);
    /// ontology.insert_term("Bar".into(), 2u32);
    ///
    /// ontology.add_parent(1u32, 2u32);
    ///
    /// assert!(ontology.hpo(2u32).unwrap().parent_ids().contains(&1u32.into()));
    /// ```
    pub fn add_parent<I: Into<HpoTermId> + Copy, J: Into<HpoTermId> + Copy>(
        &mut self,
        parent_id: I,
        child_id: J,
    ) {
        let parent = self.hpo_terms.get_unchecked_mut(parent_id.into());
        parent.add_child(child_id);

        let child = self.hpo_terms.get_unchecked_mut(child_id.into());
        child.add_parent(parent_id);
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
    pub(crate) fn add_parent_from_bytes(&mut self, bytes: &[u8]) {
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


    /// Crates and caches the `all_parents` values for every term
    ///
    /// This method can only be called once and afterwards no new terms
    /// should be added to the Ontology anymore and no new term-parent connection
    /// should be created.
    /// Since this method caches the results, rerunning it will not cause a new
    /// calculation.
    ///
    /// # Examples
    ///
    /// ```
    /// use hpo::Ontology;
    ///
    /// let mut ontology = Ontology::default();
    /// ontology.insert_term("Root".into(), 1u32);
    /// ontology.insert_term("Foo".into(), 2u32);
    /// ontology.insert_term("Bar".into(), 3u32);
    ///
    /// ontology.add_parent(1u32, 2u32);
    /// ontology.add_parent(2u32, 3u32);
    ///
    /// // At this point #3 does not have info about grandparents
    /// assert!(!ontology.hpo(3u32).unwrap().all_parent_ids().contains(&1u32.into()));
    ///
    /// ontology.create_cache();
    /// assert!(ontology.hpo(3u32).unwrap().all_parent_ids().contains(&1u32.into()));
    /// ```
    #[must_use]
    pub fn connect_all_terms(mut self) -> Builder<ConnectedTerms> {
        let term_ids: Vec<HpoTermId> = self.hpo_terms.keys();

        for id in term_ids {
            self.create_cache_of_grandparents(id);
        }
        transition_state(self)
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
        let mut res = HpoGroup::default();
        let parents = self.hpo_terms.get_unchecked(term_id).parents().clone();
        for parent in &parents {
            let grandparents = self.all_grandparents(parent);
            for gp in grandparents {
                res.insert(gp);
            }
        }
        let term = self.hpo_terms.get_unchecked_mut(term_id);
        *term.all_parents_mut() = res.bitor(&parents);
    }

    /// This method is part of the cache creation to link all terms to their
    /// direct and indirect parents (grandparents)
    ///
    /// # Panics
    ///
    /// This method will panic if the `term_id` is not present in the Ontology
    fn all_grandparents(&mut self, term_id: HpoTermId) -> &HpoGroup {
        if !self.hpo_terms.get_unchecked(term_id).parents_cached() {
            self.create_cache_of_grandparents(term_id);
        }
        let term = self.hpo_terms.get_unchecked(term_id);
        term.all_parents()
    }
}


impl Builder<ConnectedTerms> {
    /// Adds genes to the ontoloigy and connects them to connected terms
    ///
    /// This method is part of the Ontology-building, based on the binary
    /// data format and requires a specified data layout.
    ///
    /// It connects all connected terms and their parents properly. The
    /// method assumes that the bytes encode all gene-term connections.
    ///
    /// See [`Gene::as_bytes`] for explanation of the binary layout
    pub (crate) fn add_genes_from_bytes(&mut self, bytes: &[u8]) -> HpoResult<()> {
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
    pub (crate) fn add_omim_disease_from_bytes(&mut self, bytes: &[u8]) -> HpoResult<()> {
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

    /// Adds [`OrphaDisease`]s to the ontoloigy and connects them to connected terms
    ///
    /// This method is part of the Ontology-building, based on the binary
    /// data format and requires a specified data layout.
    ///
    /// It connects all connected terms and their parents properly. The
    /// method assumes that the bytes encode all Disease-term connections.
    ///
    /// See [`OrphaDisease::as_bytes`] for explanation of the binary layout
    pub (crate) fn add_orpha_disease_from_bytes(&mut self, bytes: &[u8]) -> HpoResult<()> {
        let mut idx: usize = 0;
        loop {
            if idx >= bytes.len() {
                break;
            }
            let disease_len = u32_from_bytes(&bytes[idx..]) as usize;
            let disease = OrphaDisease::try_from(&bytes[idx..idx + disease_len])?;
            for term in disease.hpo_terms() {
                self.link_orpha_disease_term(term, *disease.id())?;
            }
            self.orpha_diseases.insert(*disease.id(), disease);
            idx += disease_len;
        }
        Ok(())
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
    ///
    /// # Examples
    ///
    /// ```
    /// use hpo::Ontology;
    /// use hpo::annotations::GeneId;
    ///
    /// let mut ontology = Ontology::default();
    /// ontology.insert_term("Term-Foo".into(), 1u32);
    /// ontology.add_gene("Foo", GeneId::from(5));
    /// ontology.link_gene_term(1u32, GeneId::from(5u32)).unwrap();
    ///
    /// let term = ontology.hpo(1u32).unwrap();
    /// assert_eq!(term.genes().next().unwrap().name(), "Foo");
    /// ```
    pub fn link_gene_term<I: Into<HpoTermId>>(
        &mut self,
        term_id: I,
        gene_id: GeneId,
    ) -> HpoResult<()> {
        let term = self.hpo_terms.get_mut(term_id.into()).ok_or(HpoError::DoesNotExist)?;

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
    ///
    /// # Examples
    ///
    /// ```
    /// use hpo::Ontology;
    /// use hpo::annotations::{Disease, OmimDiseaseId};
    ///
    /// let mut ontology = Ontology::default();
    /// ontology.insert_term("Term-Foo".into(), 1u32);
    /// ontology.add_omim_disease("Foo", "5");
    /// ontology.link_omim_disease_term(1u32, OmimDiseaseId::from(5u32)).unwrap();
    ///
    /// let term = ontology.hpo(1u32).unwrap();
    /// assert_eq!(term.omim_diseases().next().unwrap().name(), "Foo");
    /// ```
    pub fn link_omim_disease_term<I: Into<HpoTermId>>(
        &mut self,
        term_id: I,
        omim_disease_id: OmimDiseaseId,
    ) -> HpoResult<()> {
        let term = self.hpo_terms.get_mut(term_id.into()).ok_or(HpoError::DoesNotExist)?;

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

    /// Add the [`OrphaDisease`] as annotation to the [`HpoTerm`]
    ///
    /// The disease will be recursively connected to all parent `HpoTerms` as well.
    ///
    /// This method does not add the HPO-term to the [`OrphaDisease`], this
    /// must be handled by the client.
    ///
    /// # Errors
    ///
    /// If the HPO term is not present, an [`HpoError`] is returned
    ///
    /// # Examples
    ///
    /// ```
    /// use hpo::Ontology;
    /// use hpo::annotations::{Disease, OrphaDiseaseId};
    ///
    /// let mut ontology = Ontology::default();
    /// ontology.insert_term("Term-Foo".into(), 1u32);
    /// ontology.add_orpha_disease("Foo", "5");
    /// ontology.link_orpha_disease_term(1u32, OrphaDiseaseId::from(5u32)).unwrap();
    ///
    /// let term = ontology.hpo(1u32).unwrap();
    /// assert_eq!(term.orpha_diseases().next().unwrap().name(), "Foo");
    /// ```
    pub fn link_orpha_disease_term<I: Into<HpoTermId>>(
        &mut self,
        term_id: I,
        orpha_disease_id: OrphaDiseaseId,
    ) -> HpoResult<()> {
        let term = self.hpo_terms.get_mut(term_id.into()).ok_or(HpoError::DoesNotExist)?;

        if term.add_orpha_disease(orpha_disease_id) {
            // If the disease is already associated to the term, this branch will
            // be skipped. That is desired, because by definition
            // all parent terms are already linked as well
            let parents = term.all_parents().clone();
            for parent in &parents {
                self.link_orpha_disease_term(parent, orpha_disease_id)?;
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
    ///
    /// # Examples
    ///
    /// ```
    /// use hpo::Ontology;
    ///
    /// let mut ontology = Ontology::default();
    ///
    /// // [all kind of logic to add terms, diseases, genes....]
    ///
    /// ontology.calculate_information_content().unwrap();
    /// ```
    #[must_use]
    pub fn calculate_information_content(mut self) -> HpoResult<Builder<FullyAnnotated>> {
        self.calculate_gene_ic()?;
        self.calculate_omim_disease_ic()?;
        self.calculate_orpha_disease_ic()?;

        Ok(transition_state(self))
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

    /// Calculates the Orpha-Disease-specific Information Content for every term
    ///
    /// If no diseases are present in the Ontology, no IC are calculated
    fn calculate_orpha_disease_ic(&mut self) -> HpoResult<()> {
        let n_orpha_diseases = self.orpha_diseases.len();

        for term in self.hpo_terms.values_mut() {
            let current_diseases = term.orpha_diseases().len();
            term.information_content_mut()
                .set_orpha_disease(n_orpha_diseases, current_diseases)?;
        }
        Ok(())
    }
}


impl Builder<FullyAnnotated> {
    pub fn build_with_defaults(self) -> HpoResult<Ontology> {
        let mut ont = Ontology { 
            hpo_terms: self.hpo_terms,
            genes: self.genes,
            omim_diseases: self.omim_diseases,
            orpha_diseases: self.orpha_diseases,
            hpo_version: self.hpo_version,
            ..Default::default()
        };
        ont.set_default_categories()?;
        ont.set_default_modifier()?;
        Ok(ont)
    }
}

impl<T> Builder<T> {
    pub fn set_hpo_version(&mut self, version: (u16, u8, u8)) {
        self.hpo_version = version;
    }

    /// Parses `Bytes` into the Jax-Ontology release version
    pub (crate) fn hpo_version_from_bytes(&mut self, bytes: &Bytes) -> HpoResult<usize> {
        if bytes.version() == BinaryVersion::V1 {
            self.set_hpo_version((0u16, 0u8, 0u8));
            Ok(0)
        } else {
            if bytes.len() < 4 {
                return Err(HpoError::ParseBinaryError);
            }
            let year = u16::from_be_bytes([bytes[0], bytes[1]]);
            let month = u8::from_be_bytes([bytes[2]]);
            let day = u8::from_be_bytes([bytes[3]]);
            self.set_hpo_version((year, month, day));
            Ok(4)
        }
    }
}


/*
struct OntologyBuilder{}
impl OntologyBuilder {
    fn add_terms(&mut self, )
}

1. add terms
2. connect terms to parents
3. add genes or diseases
4. connect genes or diseases with terms (2 must be finished)

stateDiagram-v2
    Builder --> Builder : add terms
    Builder --> Builder: add annotations (gene, disease)
    Builder --> Builder2: connect parents and children (add_parent())
    Builder2 --> Builder3: cache all terms and parents (create_cache())
    Builder2 --> Builder2: add annotations (gene, disease)
    Builder3 --> Builder3: add annotations (gene, disease)
    Builder3 --> Builder3: set_categories, set_modifier
    Builder3 --> Builder3: link annotations to terms
    Builder3 --> Builder4: calculate information content
    Builder4 --> Builder4: set_categories, set_modifier


Builder<LooseCollection>
|
add_parent()
|
V
Builder<AllTerms>
|
crate_cache()
|
V
Builder<ConnectedTerms>
|
calculate_information_content()
|
V
Builder<FullyAnnotated>
|
ontology()
|
V
Ontology
*/
