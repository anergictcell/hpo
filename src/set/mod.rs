//! An `HpoSet` can represent e.g. the clinical information of a patient or the symptoms of a disease
use crate::HpoResult;
use crate::annotations::Genes;
use crate::annotations::OmimDiseases;
use crate::similarity::GroupSimilarity;
use crate::similarity::Similarity;
use crate::similarity::SimilarityCombiner;
use crate::term::HpoGroup;
use crate::term::HpoTerms;
use crate::term::InformationContent;
use crate::HpoTerm;
use crate::Ontology;

/// A set of unique HPO terms
///
/// As in a set, each term can only appear once
/// though that is not yet guaranteed in the implementation (TODO)
pub struct HpoSet<'a> {
    ontology: &'a Ontology,
    group: HpoGroup,
}

impl<'a> HpoSet<'a> {
    /// Constructs an [`HpoSet`]
    pub fn new(ontology: &'a Ontology, group: HpoGroup) -> Self {
        Self { ontology, group }
    }

    /// Returns a new `HpoSet` that contains only the child-most terms
    ///
    /// This means that it only contains terms that don't have a child
    /// term present in the set.
    pub fn child_nodes(&mut self) -> HpoSet {
        let group = self
            .group
            .iter()
            .filter(|term1_id| {
                self.group.iter().all(|term2_id| {
                    !self
                        .ontology
                        .get_unchecked(term2_id)
                        .all_parents()
                        .contains(term1_id)
                })
            })
            .collect();
        HpoSet::new(self.ontology, group)
    }

    /// Returns the number of terms in the set
    pub fn len(&self) -> usize {
        self.group.len()
    }

    /// Returns true if there are no terms in the set
    pub fn is_empty(&self) -> bool {
        self.group.is_empty()
    }

    /// Returns a new set without modifier terms
    ///
    /// This is not yet implemented
    pub fn remove_modifier(&mut self) {
        unimplemented!()
    }

    /// Returns a new set where all obsolete terms are replaced
    ///
    /// This is not yet implemented
    pub fn replace_obsolete(&mut self, _ontology: &Ontology) {
        unimplemented!()
    }

    /// Returns all [`crate::annotations::GeneId`]s that are associated to the set
    ///
    /// # Panics
    ///
    /// When an `HpoTermId` of the set is not part of the Ontology
    pub fn gene_ids(&self) -> Genes {
        self.group
            .into_iter()
            .map(|term_id| self.ontology.get_unchecked(term_id).genes())
            .fold(Genes::default(), |acc, element| &acc | element)
    }

    /// Returns all [`crate::annotations::OmimDiseaseId`]s that are associated to the set
    ///
    /// # Panics
    ///
    /// When an `HpoTermId` of the set is not part of the Ontology
    pub fn omim_disease_ids(&self) -> OmimDiseases {
        self.group
            .into_iter()
            .map(|term_id| self.ontology.get_unchecked(term_id).omim_diseases())
            .fold(OmimDiseases::default(), |acc, element| &acc | element)
    }

    /// Calculates and returns the aggregated [`InformationContent`] of the set
    ///
    /// The `InformationContent` is not cached internally, so this operation
    /// is not cheap
    ///
    /// # Errors
    ///
    /// - When the ontology or set have more than `u16::MAX` genes or diseases
    ///
    /// # Panics
    ///
    /// - When an `HpoTermId` of the set is not part of the Ontology
    pub fn information_content(&self) -> HpoResult<InformationContent> {
        let n_diseases = self.ontology.omim_diseases().len();
        let n_genes = self.ontology.genes().len();

        let mut ic = InformationContent::default();
        ic.set_gene(n_genes, self.gene_ids().len())?;
        ic.set_omim_disease(n_diseases, self.omim_disease_ids().len())?;
        Ok(ic)
    }

    pub fn common_ancestor_ids(&self) -> HpoGroup {
        unimplemented!()
    }

    /// Returns the [`HpoTerm`] at the given index
    ///
    /// # Panics
    ///
    /// - When an `HpoTermId` of the set is not part of the Ontology
    pub fn get(&self, index: usize) -> Option<HpoTerm<'a>> {
        let term = self.group.get(index)?;
        Some(HpoTerm::try_new(self.ontology, *term).unwrap())
    }

    pub fn similarity<S: Similarity, C: SimilarityCombiner>(
        &self,
        other: &HpoSet,
        similarity: S,
        combiner: C,
    ) -> f32 {
        let group_sim = GroupSimilarity::new(combiner, similarity);
        group_sim.calculate(self, other)
    }
}

impl<'a> IntoIterator for &'a HpoSet<'a> {
    type Item = HpoTerm<'a>;
    type IntoIter = HpoTerms<'a>;
    fn into_iter(self) -> Self::IntoIter {
        HpoTerms::new(&self.group, self.ontology)
    }
}
