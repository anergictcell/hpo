//! [`Builder`] can be used to manually create custom Ontologies
//!
//! In most cases, this is not recommended, use the
//! built-in functions in [`Ontology`](`crate::Ontology`) instead.

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

use crate::HpoError;
use crate::HpoResult;
use crate::{u32_from_bytes, HpoTermId, Ontology};

use crate::ontology::termarena::Arena;

/// State of [`Builder`] that only contains some 'loose', unconnected terms
pub struct LooseCollection;

/// State of [`Builder`] that only contains all terms of the ontology,
/// but not yet connected to each other.
pub struct AllTerms;

/// State of [`Builder`] that contains all terms, connected to each other
pub struct ConnectedTerms;

/// State of [`Builder`] that contains all terms and all gene and disease annotation
pub struct FullyAnnotated;

fn transition_state<TX, TY>(builder: Builder<TX>) -> Builder<TY> {
    Builder::<TY> {
        hpo_terms: builder.hpo_terms,
        genes: builder.genes,
        omim_diseases: builder.omim_diseases,
        orpha_diseases: builder.orpha_diseases,
        hpo_version: builder.hpo_version,
        categories: builder.categories,
        modifier: builder.modifier,
        state: PhantomData,
    }
}

/// Builder to manually create an Ontology
///
/// There should rarely, if ever, be any need to build custom Ontologies.
/// The connections within the HPO, along with gene and disease annotations
/// are quite complex and it's trivially easy to mess this up, when doing it
/// manually.
///
/// To build a full ontology, the builder transitions between the following states:
///
/// ```text
/// Builder<LooseCollection> : Add individual terms to the Ontology
///   |
///   terms_complete()
///   |
///   v
/// Builder<AllTerms> : Define parent-child relationships
///   |
///   connect_all_terms()
///   |
///   v
/// Builder<ConnectedTerms>: Annotate terms and their ancestors with genes/diseases
///   |
///   calculate_information_content()
///   |
///   v
/// Builder<FullyAnnotated>
///   |
///   build_with_defaults() or build_minimal()
///   |
///   v
/// Ontology
/// ```
///
/// # Examples
///
/// ```
/// use hpo::builder::Builder;
///
/// let mut builder = Builder::new();
///
/// // Add three terms
/// builder.new_term("Root", 1u32);
/// builder.new_term("First child", 2u32);
/// builder.new_term("Second child", 3u32);
///
/// // before connecting terms, indicate that all terms have been added
/// let mut builder = builder.terms_complete();
///
/// // Connect both childs to the root term
/// builder.add_parent(1u32, 2u32);
/// builder.add_parent(1u32, 3u32);
///
/// // Build all connections and cache the connections
/// let mut builder = builder.connect_all_terms();
///
/// builder.annotate_gene(11u32.into(), "Gene1", 2u32.into()).unwrap();
/// builder.annotate_omim_disease(22u32.into(), "Disease 1", 3u32.into()).unwrap();
///
/// // Indicate that all annotations are added an calculate the information content
/// let mut builder = builder.calculate_information_content().unwrap();
///
/// // Build an Ontology
/// let ontology = builder.build_minimal();
///
/// assert_eq!(ontology.len(), 3);
///
/// let root_term = ontology.hpo(1u32).unwrap();
/// assert_eq!(root_term.name(), "Root");
/// ```
///
pub struct Builder<T> {
    hpo_terms: Arena,
    genes: HashMap<GeneId, Gene>,
    omim_diseases: HashMap<OmimDiseaseId, OmimDisease>,
    orpha_diseases: HashMap<OrphaDiseaseId, OrphaDisease>,
    hpo_version: (u16, u8, u8),
    categories: HpoGroup,
    modifier: HpoGroup,
    state: PhantomData<T>,
}

impl Default for Builder<LooseCollection> {
    fn default() -> Self {
        Self::new()
    }
}

/// In this state, the Builder contains a loose collection
/// of HPO-terms. They don't yet have any relation to each
/// other or any associated annotations
impl Builder<LooseCollection> {
    /// Creates a new `Builder` instance to manually crate an Ontology
    pub fn new() -> Builder<LooseCollection> {
        Builder::<LooseCollection> {
            hpo_terms: Arena::default(),
            genes: HashMap::default(),
            omim_diseases: HashMap::default(),
            orpha_diseases: HashMap::default(),
            hpo_version: (0u16, 0u8, 0u8),
            categories: HpoGroup::default(),
            modifier: HpoGroup::default(),
            state: PhantomData,
        }
    }

    /// Adds [`crate::HpoTerm`]s to the ontology
    ///
    /// This method is part of the Ontology-building, based on the binary
    /// data format and requires a specified data layout.
    ///
    /// The method assumes that the data is in the right format and also
    /// assumes that the caller takes care of handling all consistencies
    /// like parent-child connection etc.
    ///
    /// See [`HpoTermInternal::as_bytes`] for explanation of the binary layout.
    pub(crate) fn add_terms_from_bytes(&mut self, bytes: Bytes) {
        for term in BinaryTermBuilder::new(bytes) {
            self.add_term(term);
        }
    }

    /// Insert an `HpoTermInternal` to the ontology
    ///
    /// This method does not link the term to its parents or to any annotations.
    /// Since `HpoTermInternal` is a crate-private struct, this method
    /// is only available in-crate.
    pub(crate) fn add_term(&mut self, term: HpoTermInternal) {
        self.hpo_terms.insert(term);
    }

    /// Adds a new term to the ontology
    ///
    /// The term does not have any connections to other terms or any gene/disease
    /// annotations. Parents and children of terms can be added in the
    /// `Builder<AllTerms>` state.
    ///
    /// # Examples
    ///
    /// ```
    /// use hpo::builder::Builder;
    ///
    /// let mut builder = Builder::new();
    ///
    /// // Add three terms
    /// builder.new_term("Root", 1u32);
    /// builder.new_term("First child", 2u32);
    /// builder.new_term("Second child", 3u32);
    ///
    /// // quickly transition through all stages to build the ontology
    /// let ontology = builder
    ///     .terms_complete()
    ///     .connect_all_terms()
    ///     .calculate_information_content().unwrap()
    ///     .build_minimal();
    ///
    /// assert_eq!(ontology.len(), 3);
    /// ```
    pub fn new_term<I: Into<HpoTermId>>(&mut self, name: &str, id: I) {
        let term = HpoTermInternal::new(name.to_string(), id.into());
        self.add_term(term);
    }

    /// Transitions the state to `Builder<AllTerms>`
    ///
    /// This method indicates that all terms have been added. It is not possible
    /// to add new terms afterwards.
    /// Transitioning to `Builder<AllTerms>` is required to crate parent-child
    /// connections between terms.
    #[must_use]
    pub fn terms_complete(self) -> Builder<AllTerms> {
        transition_state(self)
    }
}

impl Builder<AllTerms> {
    /// Add a connection from an [`HpoTerm`](`crate::HpoTerm`) to its parent
    ///
    /// This method is called once for every dependency in the Ontology during the initialization.
    ///
    /// # Errors
    ///
    /// This method will return `HpoError::DoesNotExist` when the parent or child
    /// `HpoTerm` does not exist.
    ///
    /// # Examples
    ///
    /// ```
    /// use hpo::builder::Builder;
    /// # use hpo::builder::AllTerms;
    ///
    /// fn example_builder() -> Builder<AllTerms>
    /// # {
    /// # let mut builder = Builder::new();
    /// # builder.new_term("Foo", 1u32);
    /// # builder.new_term("Bar", 2u32);
    /// # builder.terms_complete()
    /// # }
    ///
    /// let mut builder: Builder<AllTerms> = example_builder();
    ///
    /// // connect a term to its parent
    /// builder.add_parent(1u32, 2u32).unwrap();
    ///
    /// // quickly transition through all stages to build the ontology
    /// let ontology = builder
    ///     .connect_all_terms()
    ///     .calculate_information_content().unwrap()
    ///     .build_minimal();
    ///
    /// assert!(ontology.hpo(2u32).unwrap().parent_ids().contains(&1u32.into()));
    /// ```
    pub fn add_parent<I: Into<HpoTermId> + Copy, J: Into<HpoTermId> + Copy>(
        &mut self,
        parent_id: I,
        child_id: J,
    ) -> HpoResult<()> {
        let parent = self
            .hpo_terms
            .get_mut(parent_id.into())
            .ok_or(HpoError::DoesNotExist)?;
        parent.add_child(child_id);

        let child = self
            .hpo_terms
            .get_mut(child_id.into())
            .ok_or(HpoError::DoesNotExist)?;
        child.add_parent(parent_id);
        Ok(())
    }

    /// Add a connection from an [`HpoTerm`](`crate::HpoTerm`) to its parent
    ///
    /// This method is called once for every dependency in the Ontology during the initialization.
    ///
    /// # Panics
    ///
    /// This method will panic if the `parent_id` or `child_id` is not present in the Ontology
    ///
    pub(crate) fn add_parent_unchecked<I: Into<HpoTermId> + Copy, J: Into<HpoTermId> + Copy>(
        &mut self,
        parent_id: I,
        child_id: J,
    ) {
        let parent = self.hpo_terms.get_unchecked_mut(parent_id.into());
        parent.add_child(child_id);

        let child = self.hpo_terms.get_unchecked_mut(child_id.into());
        child.add_parent(parent_id);
    }

    /// Connects an [`HpoTerm`](`crate::HpoTerm`) to its parent term
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
    /// to the contained data or if a `parent_id` or `child_id` is not present
    /// in the Ontology
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
                self.add_parent_unchecked(parent, term);
                idx += 4;
            }
        }
    }

    /// Transitions the state to `Builder<ConnectedTerms>`
    ///
    /// After changing the state no new terms can be added to the Ontology
    /// anymore and no new term-parent connection can should be created.
    ///
    /// # Examples
    ///
    /// ```
    /// use hpo::builder::Builder;
    /// # use hpo::builder::{AllTerms, ConnectedTerms};
    ///
    /// fn example_builder() -> Builder<AllTerms>
    /// # {
    /// # let mut builder = Builder::new();
    /// # builder.new_term("Foo", 1u32);
    /// # builder.new_term("Bar", 2u32);
    /// # let mut builder = builder.terms_complete();
    /// # builder.add_parent(1u32, 2u32).unwrap();
    /// # builder
    /// # }
    ///
    /// let mut builder: Builder<AllTerms> = example_builder();
    ///
    /// // connect all the terms and return a `Builder<ConnectedTerms>`
    /// let builder: Builder<ConnectedTerms> = builder.connect_all_terms();
    ///
    /// // quickly transition through all stages to build the ontology
    /// let ontology = builder
    ///     .calculate_information_content().unwrap()
    ///     .build_minimal();
    ///
    /// assert!(ontology.hpo(2u32).unwrap().parent_ids().contains(&1u32.into()));
    /// ```
    #[must_use]
    pub fn connect_all_terms(mut self) -> Builder<ConnectedTerms> {
        for id in self.hpo_terms.keys() {
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
    /// Add the [`Gene`] as annotation to the [`HpoTerm`](`crate::HpoTerm`)
    ///
    /// The gene will be recursively connected to all parent `HpoTerms` as well.
    ///
    /// # Errors
    ///
    /// If the HPO term is not present, an [`HpoError::DoesNotExist`] is returned
    ///
    /// # Examples
    ///
    /// ```
    /// use hpo::HpoTermId;
    /// use hpo::annotations::GeneId;
    /// use hpo::builder::Builder;
    /// # use hpo::builder::ConnectedTerms;
    ///
    /// fn example_builder() -> Builder<ConnectedTerms>
    /// # {
    /// # let mut builder = Builder::new();
    /// # builder.new_term("Foo", 1u32);
    /// # builder.new_term("Bar", 2u32);
    /// # let mut builder = builder.terms_complete();
    /// # builder.add_parent(1u32, 2u32).unwrap();
    /// # builder.connect_all_terms()
    /// # }
    ///
    /// let mut builder: Builder<ConnectedTerms> = example_builder();
    ///
    /// builder.annotate_gene(GeneId::from(5), "Gene 1", HpoTermId::from(1u32));
    ///
    /// // quickly transition through all stages to build the ontology
    /// let ontology = builder
    ///     .calculate_information_content().unwrap()
    ///     .build_minimal();
    ///    
    /// let term = ontology.hpo(1u32).unwrap();
    /// assert!(term.genes().find(|gene| gene.name() == "Gene 1").is_some());
    /// assert!(term.genes().find(|gene| gene.name() == "Foobar").is_none());
    /// ```
    #[allow(clippy::missing_panics_doc)]
    pub fn annotate_gene(
        &mut self,
        gene_id: GeneId,
        gene_name: &str,
        term_id: HpoTermId,
    ) -> HpoResult<()> {
        self.add_gene(gene_name, gene_id);
        let gene = self
            .genes
            .get_mut(&gene_id)
            .expect("Gene is present because it was just add_omim_disease");

        gene.add_term(term_id);
        self.link_gene_term(term_id, gene_id)?;
        Ok(())
    }

    /// Add the [`OmimDisease`] as annotation to the [`HpoTerm`](`crate::HpoTerm`)
    ///
    /// The disease will be recursively connected to all parent `HpoTerms` as well.
    ///
    /// # Errors
    ///
    /// If the HPO term is not present, an [`HpoError::DoesNotExist`] is returned
    ///
    /// # Examples
    ///
    /// ```
    /// use hpo::HpoTermId;
    /// use hpo::annotations::{Disease, OmimDiseaseId};
    /// use hpo::builder::Builder;
    /// # use hpo::builder::ConnectedTerms;
    ///
    /// fn example_builder() -> Builder<ConnectedTerms>
    /// # {
    /// # let mut builder = Builder::new();
    /// # builder.new_term("Foo", 1u32);
    /// # builder.new_term("Bar", 2u32);
    /// # let mut builder = builder.terms_complete();
    /// # builder.add_parent(1u32, 2u32).unwrap();
    /// # builder.connect_all_terms()
    /// # }
    ///
    /// let mut builder: Builder<ConnectedTerms> = example_builder();
    ///
    /// builder.annotate_omim_disease(OmimDiseaseId::from(5), "Disease 1", HpoTermId::from(1u32));
    ///
    /// // quickly transition through all stages to build the ontology
    /// let ontology = builder
    ///     .calculate_information_content().unwrap()
    ///     .build_minimal();
    ///    
    /// let term = ontology.hpo(1u32).unwrap();
    /// assert!(term.omim_diseases().find(|disease| disease.name() == "Disease 1").is_some());
    /// assert!(term.omim_diseases().find(|disease| disease.name() == "Foobar").is_none());
    /// ```
    #[allow(clippy::missing_panics_doc)]
    pub fn annotate_omim_disease(
        &mut self,
        omim_id: OmimDiseaseId,
        omim_name: &str,
        term_id: HpoTermId,
    ) -> HpoResult<()> {
        self.add_omim_disease(omim_name, omim_id);
        let gene = self
            .omim_diseases
            .get_mut(&omim_id)
            .expect("Gene is present because it was just add_omim_disease");

        gene.add_term(term_id);
        self.link_omim_disease_term(term_id, omim_id)?;

        Ok(())
    }

    /// Add the [`OrphaDisease`] as annotation to the [`HpoTerm`](`crate::HpoTerm`)
    ///
    /// The disease will be recursively connected to all parent `HpoTerms` as well.
    ///
    /// # Errors
    ///
    /// If the HPO term is not present, an [`HpoError::DoesNotExist`] is returned
    ///
    /// # Examples
    ///
    /// ```
    /// use hpo::HpoTermId;
    /// use hpo::annotations::{Disease, OrphaDiseaseId};
    /// use hpo::builder::Builder;
    /// # use hpo::builder::ConnectedTerms;
    ///
    /// fn example_builder() -> Builder<ConnectedTerms>
    /// # {
    /// # let mut builder = Builder::new();
    /// # builder.new_term("Foo", 1u32);
    /// # builder.new_term("Bar", 2u32);
    /// # let mut builder = builder.terms_complete();
    /// # builder.add_parent(1u32, 2u32).unwrap();
    /// # builder.connect_all_terms()
    /// # }
    ///
    /// let mut builder: Builder<ConnectedTerms> = example_builder();
    ///
    /// builder.annotate_orpha_disease(OrphaDiseaseId::from(5), "Disease 1", HpoTermId::from(1u32));
    ///
    /// // quickly transition through all stages to build the ontology
    /// let ontology = builder
    ///     .calculate_information_content().unwrap()
    ///     .build_minimal();
    ///    
    /// let term = ontology.hpo(1u32).unwrap();
    /// assert!(term.orpha_diseases().find(|disease| disease.name() == "Disease 1").is_some());
    /// assert!(term.orpha_diseases().find(|disease| disease.name() == "Foobar").is_none());
    /// ```
    #[allow(clippy::missing_panics_doc)]
    pub fn annotate_orpha_disease(
        &mut self,
        orpha_id: OrphaDiseaseId,
        orpha_name: &str,
        term_id: HpoTermId,
    ) -> HpoResult<()> {
        self.add_orpha_disease(orpha_name, orpha_id);
        let gene = self
            .orpha_diseases
            .get_mut(&orpha_id)
            .expect("Gene is present because it was just add_orpha_disease");

        gene.add_term(term_id);
        self.link_orpha_disease_term(term_id, orpha_id)?;

        Ok(())
    }

    /// Calculates the [`crate::term::InformationContent`]s for every term
    /// and transitions to the `FullyAnnotated` state
    ///
    /// This method should only be called **after** all terms are added,
    /// connected and all genes and diseases are linked as well.
    ///
    /// # Errors
    ///
    /// This method returns an error if there are more Genes or Terms than `u16::MAX`
    /// because larger numbers can't be safely converted to `f32`
    ///
    /// # Examples
    ///
    /// ```
    /// use hpo::HpoTermId;
    /// use hpo::annotations::GeneId;
    /// use hpo::builder::Builder;
    /// # use hpo::builder::ConnectedTerms;
    ///
    /// fn example_builder() -> Builder<ConnectedTerms>
    /// # {
    /// # let mut builder = Builder::new();
    /// # builder.new_term("Foo", 1u32);
    /// # builder.new_term("Bar", 2u32);
    /// # let mut builder = builder.terms_complete();
    /// # builder.add_parent(1u32, 2u32).unwrap();
    /// # let mut builder = builder.connect_all_terms();
    /// # builder
    /// # }
    ///
    /// let mut builder: Builder<ConnectedTerms> = example_builder();
    /// builder.annotate_gene(GeneId::from(1), "Gene 1", HpoTermId::from(1u32));
    /// builder.annotate_gene(GeneId::from(2), "Gene 2", HpoTermId::from(2u32));
    ///
    /// // transition to final state of `Builder`
    /// let builder = builder.calculate_information_content().unwrap();
    ///
    /// let ontology = builder.build_minimal();
    ///
    /// let gene_ic = ontology
    ///     .hpo(2u32).unwrap()
    ///     .information_content()
    ///     .gene();
    /// assert!(gene_ic > 0.0, "{gene_ic}");
    /// ```
    pub fn calculate_information_content(mut self) -> HpoResult<Builder<FullyAnnotated>> {
        self.calculate_gene_ic()?;
        self.calculate_omim_disease_ic()?;
        self.calculate_orpha_disease_ic()?;

        Ok(transition_state(self))
    }

    /// Adds a [`Gene`](`crate::annotations::Gene`) to the ontology
    ///
    /// The gene is not yet linked to any terms, this must be done
    /// through [`Builder<ConnectedTerms>::annotate_gene`](`Builder::annotate_gene`)
    ///
    /// # Note:
    ///
    /// There is rarely need to call this method directly. The preferred way is to use
    /// [`Builder<ConnectedTerms>::annotate_gene`](`Builder::annotate_gene`)
    pub fn add_gene(&mut self, gene_name: &str, gene_id: GeneId) {
        if let Entry::Vacant(entry) = self.genes.entry(gene_id) {
            entry.insert(Gene::new(gene_id, gene_name));
        }
    }

    /// Adds an [`OmimDisease`](`crate::annotations::OmimDisease`) to the ontology
    ///
    /// The gene is not yet linked to any terms, this must be done
    /// through [`Builder<ConnectedTerms>::annotate_omim_disease`](`Builder::annotate_omim_disease`)
    ///
    /// # Note:
    ///
    /// There is rarely need to call this method directly. The preferred way is to use
    /// [`Builder<ConnectedTerms>::annotate_omim_disease`](`Builder::annotate_omim_disease`)
    pub fn add_omim_disease(
        &mut self,
        omim_disease_name: &str,
        omim_disease_id: OmimDiseaseId,
    ) -> OmimDiseaseId {
        match self.omim_diseases.entry(omim_disease_id) {
            std::collections::hash_map::Entry::Occupied(_) => omim_disease_id,
            std::collections::hash_map::Entry::Vacant(entry) => {
                entry.insert(OmimDisease::new(omim_disease_id, omim_disease_name));
                omim_disease_id
            }
        }
    }

    /// Adds an [`OrphaDisease`](`crate::annotations::OrphaDisease`) to the ontology
    ///
    /// The gene is not yet linked to any terms, this must be done
    /// through [`Builder<ConnectedTerms>::annotate_orpha_disease`](`Builder::annotate_orpha_disease`)
    ///
    /// # Note:
    ///
    /// There is rarely need to call this method directly. The preferred way is to use
    /// [`Builder<ConnectedTerms>::annotate_orpha_disease`](`Builder::annotate_orpha_disease`)
    pub fn add_orpha_disease(
        &mut self,
        orpha_disease_name: &str,
        orpha_disease_id: OrphaDiseaseId,
    ) -> OrphaDiseaseId {
        match self.orpha_diseases.entry(orpha_disease_id) {
            std::collections::hash_map::Entry::Occupied(_) => orpha_disease_id,
            std::collections::hash_map::Entry::Vacant(entry) => {
                entry.insert(OrphaDisease::new(orpha_disease_id, orpha_disease_name));
                orpha_disease_id
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
    pub(crate) fn add_genes_from_bytes(&mut self, bytes: &[u8]) -> HpoResult<()> {
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
    pub(crate) fn add_omim_disease_from_bytes(&mut self, bytes: &[u8]) -> HpoResult<()> {
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
    pub(crate) fn add_orpha_disease_from_bytes(&mut self, bytes: &[u8]) -> HpoResult<()> {
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

    /// Add the [`Gene`] as annotation to the [`HpoTerm`](`crate::HpoTerm`)
    ///
    /// The gene will be recursively connected to all parent `HpoTerms` as well.
    ///
    /// This method does not add the HPO-term to the [`Gene`], this
    /// must be handled by the client.
    ///
    /// # Errors
    ///
    /// If the HPO term is not present, an [`HpoError::DoesNotExist`] is returned
    ///
    fn link_gene_term(&mut self, term_id: HpoTermId, gene_id: GeneId) -> HpoResult<()> {
        let term = self
            .hpo_terms
            .get_mut(term_id)
            .ok_or(HpoError::DoesNotExist)?;

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

    /// Add the [`OmimDisease`] as annotation to the [`HpoTerm`](`crate::HpoTerm`)
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
    fn link_omim_disease_term(
        &mut self,
        term_id: HpoTermId,
        omim_disease_id: OmimDiseaseId,
    ) -> HpoResult<()> {
        let term = self
            .hpo_terms
            .get_mut(term_id)
            .ok_or(HpoError::DoesNotExist)?;

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

    /// Add the [`OrphaDisease`] as annotation to the [`HpoTerm`](`crate::HpoTerm`)
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
    fn link_orpha_disease_term(
        &mut self,
        term_id: HpoTermId,
        orpha_disease_id: OrphaDiseaseId,
    ) -> HpoResult<()> {
        let term = self
            .hpo_terms
            .get_mut(term_id)
            .ok_or(HpoError::DoesNotExist)?;

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
    /// Builds the [`Ontology`] with default settings
    ///
    /// This method can only be used with the standard HPO from Jax
    /// and will most likely not work with custom ontologies.
    ///
    /// # Errors
    ///
    /// This method requires that the main-category terms:
    ///
    /// - `HP:0000001 | All`
    /// - `HP:0000118 | Phenotypic abnormality`
    ///
    /// are present in the Ontology.
    pub fn build_with_defaults(self) -> HpoResult<Ontology> {
        let mut ont = self.build_minimal();
        ont.set_default_categories()?;
        ont.set_default_modifier()?;
        Ok(ont)
    }

    /// Builds the [`Ontology`]
    ///
    /// This method will not specify different phenotype
    /// categories or modifier terms.
    ///
    /// Use this method only with custom ontologies. When using the standard
    /// Jax ontology, use the recommended [`Builder::build_with_defaults`]
    /// method.
    pub fn build_minimal(self) -> Ontology {
        Ontology {
            hpo_terms: self.hpo_terms,
            genes: self.genes,
            omim_diseases: self.omim_diseases,
            orpha_diseases: self.orpha_diseases,
            hpo_version: self.hpo_version,
            ..Default::default()
        }
    }
}

impl<T> Builder<T> {
    /// Defines the HPO version of the Ontology
    /// The version should be specified as \[YEAR\]-\[MONTH\]-\[DAY\], e.g.
    /// `2024-08-21`
    pub fn set_hpo_version(&mut self, version: (u16, u8, u8)) {
        self.hpo_version = version;
    }

    /// Parses `Bytes` into the Jax-Ontology release version
    ///
    /// # Errors
    ///
    /// - Wrong HPO-version format: [`HpoError::ParseBinaryError`]
    pub(crate) fn hpo_version_from_bytes(&mut self, bytes: &Bytes) -> HpoResult<usize> {
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
