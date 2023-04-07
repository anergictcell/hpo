//! An `HpoSet` can represent e.g. the clinical information of a patient or the symptoms of a disease
use crate::annotations::Genes;
use crate::annotations::OmimDiseases;
use crate::similarity::GroupSimilarity;
use crate::similarity::Similarity;
use crate::similarity::SimilarityCombiner;
use crate::term::HpoGroup;
use crate::HpoTermId;

use crate::term::{InformationContent, Iter};
use crate::HpoResult;
use crate::HpoTerm;

use crate::Ontology;

/// A set of unique HPO terms
///
/// It provides several convinience functions to operate on a set of terms.
/// A typical use-case for an [`HpoSet`] is to record the clinical information
/// of a patient. You can compare the aggregated information of the patient
/// with other patients, genes or dieases.
///
/// As in a set, each term can only appear once.
///
/// # Examples
///
/// ```
/// use hpo::term::InformationContentKind;
/// use hpo::{Ontology, HpoSet};
/// use hpo::term::HpoGroup;
/// use hpo::similarity::{Builtins, StandardCombiner};
///
/// let ontology = Ontology::from_binary("tests/example.hpo").unwrap();
///
/// // create one set
/// let mut hpos = HpoGroup::new();
/// hpos.insert(707u32);
/// hpos.insert(12639u32);
/// hpos.insert(12638u32);
/// hpos.insert(818u32);
/// hpos.insert(2715u32);
/// let set = HpoSet::new(&ontology, hpos);
/// assert_eq!(set.len(), 5);
///
/// // create another set
/// let mut hpos_2 = HpoGroup::new();
/// hpos_2.insert(100547u32);
/// hpos_2.insert(12638u32);
/// hpos_2.insert(864u32);
/// hpos_2.insert(25454u32);
/// let set_2 = HpoSet::new(&ontology, hpos_2);
/// assert_eq!(set_2.len(), 4);
///
/// let similarity = set.similarity(
///     &set_2,
///     Builtins::GraphIc(InformationContentKind::Omim),
///     StandardCombiner::default()
/// );
///
/// assert_eq!(similarity, 0.8177036);
/// ```
#[must_use]
pub struct HpoSet<'a> {
    ontology: &'a Ontology,
    group: HpoGroup,
}

impl<'a> HpoSet<'a> {
    /// Constructs an [`HpoSet`]
    ///
    /// # Examples
    ///
    /// ```
    /// use hpo::{Ontology, HpoSet};
    /// use hpo::term::HpoGroup;
    ///
    /// let ontology = Ontology::from_binary("tests/example.hpo").unwrap();
    ///
    /// let mut hpos = HpoGroup::new();
    /// hpos.insert(707u32);
    /// hpos.insert(12639u32);
    /// hpos.insert(12638u32);
    /// hpos.insert(818u32);
    /// hpos.insert(2715u32);
    /// let set = HpoSet::new(&ontology, hpos);
    /// assert_eq!(set.len(), 5);
    /// ```
    pub fn new(ontology: &'a Ontology, group: HpoGroup) -> Self {
        Self { ontology, group }
    }

    /// Returns a new `HpoSet` that contains only the child-most terms
    ///
    /// This means that it only contains terms that don't have a child
    /// term present in the set.
    ///
    /// # Examples
    ///
    /// ```
    /// use hpo::{Ontology, HpoSet};
    /// use hpo::term::HpoGroup;
    ///
    /// let ontology = Ontology::from_binary("tests/example.hpo").unwrap();
    ///
    /// let mut hpos = HpoGroup::new();
    /// hpos.insert(707u32);
    /// hpos.insert(12639u32);
    /// hpos.insert(12638u32);
    /// hpos.insert(818u32);
    /// hpos.insert(2715u32);
    /// let mut set = HpoSet::new(&ontology, hpos);
    /// assert_eq!(set.len(), 5);
    /// let children = set.child_nodes();
    /// assert_eq!(children.len(), 4);
    /// ```
    pub fn child_nodes(&mut self) -> HpoSet {
        let group = self
            .group
            .iter()
            .filter(|term1_id| {
                self.group.iter().all(|term2_id| {
                    !self
                        .ontology
                        .get(term2_id)
                        .expect("HpoTermId must be in Ontology")
                        .all_parents()
                        .contains(term1_id)
                })
            })
            .collect();
        HpoSet::new(self.ontology, group)
    }

    /// Returns the number of terms in the set
    ///
    /// # Examples
    ///
    /// ```
    /// use hpo::{Ontology, HpoSet};
    /// use hpo::term::HpoGroup;
    ///
    /// let ontology = Ontology::from_binary("tests/example.hpo").unwrap();
    ///
    /// let mut hpos = HpoGroup::new();
    /// hpos.insert(707u32);
    /// hpos.insert(12639u32);
    /// hpos.insert(12638u32);
    /// hpos.insert(818u32);
    /// hpos.insert(2715u32);
    /// let set = HpoSet::new(&ontology, hpos);
    /// assert_eq!(set.len(), 5);
    /// ```
    pub fn len(&self) -> usize {
        self.group.len()
    }

    /// Returns true if there are no terms in the set
    ///
    /// # Examples
    ///
    /// ```
    /// use hpo::{Ontology, HpoSet};
    /// use hpo::term::HpoGroup;
    ///
    /// let ontology = Ontology::from_binary("tests/example.hpo").unwrap();
    ///
    /// let hpos = HpoGroup::new();
    /// let set = HpoSet::new(&ontology, hpos);
    /// assert!(set.is_empty());
    /// ```
    pub fn is_empty(&self) -> bool {
        self.group.is_empty()
    }

    /// Removes all modifier terms in-place
    ///
    ///
    /// # Examples
    ///
    /// ```
    /// use hpo::{Ontology, HpoSet};
    /// # use hpo::term::HpoGroup;
    /// # fn method_that_returns_an_hposet<'a>(ontology: &'a Ontology) -> HpoSet<'a> {
    /// # let mut hpos = HpoGroup::new();
    /// # hpos.insert(707u32.into());
    /// # hpos.insert(12639u32.into());
    /// # hpos.insert(12638u32.into());
    /// # hpos.insert(3581u32.into());
    /// # hpos.insert(7u32.into());
    /// # HpoSet::new(ontology, hpos)
    /// # }
    ///
    /// let ontology = Ontology::from_binary("tests/example.hpo").unwrap();
    ///
    /// let mut set: HpoSet = method_that_returns_an_hposet(&ontology);
    /// assert_eq!(set.len(), 5);
    ///
    /// set.remove_modifier();
    /// assert_eq!(set.len(), 3);
    /// ```
    pub fn remove_modifier(&mut self) {
        let group: HpoGroup = self.iter().filter(|term| !term.is_modifier()).collect();
        self.group = group;
    }

    /// Returns a new set without modifier terms
    ///
    /// # Examples
    ///
    /// ```
    /// use hpo::{Ontology, HpoSet};
    /// # use hpo::term::HpoGroup;
    /// # fn method_that_returns_an_hposet<'a>(ontology: &'a Ontology) -> HpoSet<'a> {
    /// # let mut hpos = HpoGroup::new();
    /// # hpos.insert(707u32.into());
    /// # hpos.insert(12639u32.into());
    /// # hpos.insert(12638u32.into());
    /// # hpos.insert(3581u32.into());
    /// # hpos.insert(7u32.into());
    /// # HpoSet::new(ontology, hpos)
    /// # }
    ///
    /// let ontology = Ontology::from_binary("tests/example.hpo").unwrap();
    ///
    /// let set_1: HpoSet = method_that_returns_an_hposet(&ontology);
    /// assert_eq!(set_1.len(), 5);
    ///
    /// let set_2 = set_1.without_modifier();
    /// assert_eq!(set_2.len(), 3);
    ///
    /// for term in set_1.iter() {
    ///     if !set_2.contains(&term.id()) {
    ///         println!("Modifier: {} | {}", term.id(), term.name());
    ///     }
    /// }
    /// // "HP:0000007 | Autosomal recessive inheritance"
    /// // "HP:0003581 | Adult onset"
    /// ```
    pub fn without_modifier(&self) -> Self {
        let group: HpoGroup = self.iter().filter(|term| !term.is_modifier()).collect();
        Self {
            ontology: self.ontology,
            group,
        }
    }

    /// Removes all obsolete terms in-place
    ///
    /// # Panics
    ///
    /// When an `HpoTermId` of the set is not part of the Ontology
    ///
    /// # Examples
    ///
    /// The example Ontology does not contain obsolete terms,
    /// so this example does not show an actual effect.
    ///
    /// ```
    /// use hpo::{Ontology, HpoSet};
    /// use hpo::term::HpoGroup;
    /// # fn method_that_returns_an_hposet<'a>(ontology: &'a Ontology) -> HpoSet<'a> {
    /// # let mut hpos = HpoGroup::new();
    /// # hpos.insert(707u32);
    /// # hpos.insert(12639u32);
    /// # hpos.insert(12638u32);
    /// # hpos.insert(818u32);
    /// # hpos.insert(2715u32);
    /// # HpoSet::new(ontology, hpos)
    /// # }
    ///
    /// let ontology = Ontology::from_binary("tests/example.hpo").unwrap();
    ///
    /// let mut set: HpoSet = method_that_returns_an_hposet(&ontology);
    /// assert_eq!(set.len(), 5);
    ///
    /// set.remove_obsolete();
    /// assert_eq!(set.len(), 5);
    /// ```
    pub fn remove_obsolete(&mut self) {
        let group: HpoGroup = self
            .group
            .iter()
            .filter(|term_id| {
                !self
                    .ontology
                    .get(*term_id)
                    .expect("HpoTermId must be in Ontology")
                    .obsolete()
            })
            .collect();
        self.group = group;
    }

    /// Returns a new set without obsolete terms
    ///
    /// # Panics
    ///
    /// When an `HpoTermId` of the set is not part of the Ontology
    ///
    /// # Examples
    ///
    /// The example Ontology does not contain obsolete terms,
    /// so this example does not show an actual effect.
    ///
    /// ```
    /// use hpo::{Ontology, HpoSet};
    /// use hpo::term::HpoGroup;
    /// # fn method_that_returns_an_hposet<'a>(ontology: &'a Ontology) -> HpoSet<'a> {
    /// # let mut hpos = HpoGroup::new();
    /// # hpos.insert(707u32);
    /// # hpos.insert(12639u32);
    /// # hpos.insert(12638u32);
    /// # hpos.insert(818u32);
    /// # hpos.insert(2715u32);
    /// # HpoSet::new(ontology, hpos)
    /// # }
    ///
    /// let ontology = Ontology::from_binary("tests/example.hpo").unwrap();
    ///
    /// let set: HpoSet = method_that_returns_an_hposet(&ontology);
    /// assert_eq!(set.len(), 5);
    ///
    /// let new_set = set.without_obsolete();
    /// assert_eq!(new_set.len(), 5);
    /// ```
    pub fn without_obsolete(&self) -> Self {
        let group: HpoGroup = self
            .group
            .iter()
            .filter(|term_id| {
                !self
                    .ontology
                    .get(*term_id)
                    .expect("HpoTermId must be in Ontology")
                    .obsolete()
            })
            .collect();
        Self {
            ontology: self.ontology,
            group,
        }
    }

    /// Returns a new set where obsolete terms are replaced
    ///
    /// All terms with a `replaced_by` attribute are replaced accordingly.
    ///
    /// See [`HpoTerm::replaced_by`] for more information
    ///
    /// # Panics
    ///
    /// When an `HpoTermId` of the set is not part of the Ontology
    ///
    /// # Examples
    ///
    /// The example Ontology does not contain replacement terms,
    /// so this example does not show an actual effect.
    ///
    /// ```
    /// use hpo::{Ontology, HpoSet};
    /// use hpo::term::HpoGroup;
    /// # fn method_that_returns_an_hposet<'a>(ontology: &'a Ontology) -> HpoSet<'a> {
    /// # let mut hpos = HpoGroup::new();
    /// # hpos.insert(707u32);
    /// # hpos.insert(12639u32);
    /// # hpos.insert(12638u32);
    /// # hpos.insert(818u32);
    /// # hpos.insert(2715u32);
    /// # HpoSet::new(ontology, hpos)
    /// # }
    ///
    /// let ontology = Ontology::from_binary("tests/example.hpo").unwrap();
    ///
    /// let set: HpoSet = method_that_returns_an_hposet(&ontology);
    /// assert_eq!(set.len(), 5);
    ///
    /// let new_set = set.with_replaced_obsolete();
    /// assert_eq!(new_set.len(), 5);
    /// ```
    pub fn with_replaced_obsolete(&self) -> Self {
        let group: HpoGroup = self
            .group
            .iter()
            .map(|term_id| {
                self.ontology
                    .get(term_id)
                    .expect("HpoTermId must be in Ontology")
                    .replacement()
                    .unwrap_or(term_id)
            })
            .collect();
        Self {
            ontology: self.ontology,
            group,
        }
    }

    /// Replaces obsolete terms in-place
    ///
    /// All terms with a `replaced_by` attribute are replaced accordingly.
    ///
    /// See [`HpoTerm::replaced_by`] for more information
    ///
    /// # Panics
    ///
    /// When an `HpoTermId` of the set is not part of the Ontology
    ///
    /// # Examples
    ///
    /// The example Ontology does not contain replacement terms,
    /// so this example does not show an actual effect.
    ///
    /// ```
    /// use hpo::{Ontology, HpoSet};
    /// use hpo::term::HpoGroup;
    /// # fn method_that_returns_an_hposet<'a>(ontology: &'a Ontology) -> HpoSet<'a> {
    /// # let mut hpos = HpoGroup::new();
    /// # hpos.insert(707u32);
    /// # hpos.insert(12639u32);
    /// # hpos.insert(12638u32);
    /// # hpos.insert(818u32);
    /// # hpos.insert(2715u32);
    /// # HpoSet::new(ontology, hpos)
    /// # }
    ///
    /// let ontology = Ontology::from_binary("tests/example.hpo").unwrap();
    ///
    /// let mut set: HpoSet = method_that_returns_an_hposet(&ontology);
    /// assert_eq!(set.len(), 5);
    ///
    /// set.replace_obsolete();
    /// assert_eq!(set.len(), 5);
    /// ```
    pub fn replace_obsolete(&mut self) {
        let group: HpoGroup = self
            .group
            .iter()
            .map(|term_id| {
                self.ontology
                    .get(term_id)
                    .expect("HpoTermId must be in Ontology")
                    .replacement()
                    .unwrap_or(term_id)
            })
            .collect();
        self.group = group;
    }

    /// Returns all [`crate::annotations::GeneId`]s that are associated to the set
    ///
    /// # Panics
    ///
    /// When an `HpoTermId` of the set is not part of the Ontology
    ///
    /// # Examples
    ///
    /// ```
    /// use hpo::{Ontology, HpoSet};
    /// use hpo::term::HpoGroup;
    ///
    /// let ontology = Ontology::from_binary("tests/example.hpo").unwrap();
    ///
    /// let mut hpos = HpoGroup::new();
    /// hpos.insert(12639u32);
    /// hpos.insert(818u32);
    /// let set = HpoSet::new(&ontology, hpos);
    /// // `ABCD1 (HGNC:215)` is linked to `HP:0012639`
    /// assert!(set.gene_ids().contains(&215u32.into()));
    /// ```
    pub fn gene_ids(&self) -> Genes {
        self.group
            .iter()
            .map(|term_id| {
                self.ontology
                    .get(term_id)
                    .expect("HpoTermId must be in Ontology")
                    .genes()
            })
            .fold(Genes::default(), |acc, element| &acc | element)
    }

    /// Returns all [`crate::annotations::OmimDiseaseId`]s that are associated to the set
    ///
    /// # Panics
    ///
    /// When an `HpoTermId` of the set is not part of the Ontology
    ///
    /// # Examples
    ///
    /// ```
    /// use hpo::{Ontology, HpoSet};
    /// use hpo::term::HpoGroup;
    ///
    /// let ontology = Ontology::from_binary("tests/example.hpo").unwrap();
    ///
    /// let mut hpos = HpoGroup::new();
    /// hpos.insert(12639u32);
    /// hpos.insert(818u32);
    /// let set = HpoSet::new(&ontology, hpos);
    /// // `Microphthalmia, syndromic 6 (OMIM:607932)` is linked to `HP:0012639`
    /// assert!(set.omim_disease_ids().contains(&607932u32.into()));
    /// ```
    pub fn omim_disease_ids(&self) -> OmimDiseases {
        self.group
            .iter()
            .map(|term_id| {
                self.ontology
                    .get(term_id)
                    .expect("HpoTermId must be in Ontology")
                    .omim_diseases()
            })
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
    ///
    /// # Examples
    ///
    /// ```
    /// use hpo::{Ontology, HpoSet};
    /// use hpo::term::HpoGroup;
    ///
    /// let ontology = Ontology::from_binary("tests/example.hpo").unwrap();
    ///
    /// let mut hpos = HpoGroup::new();
    /// hpos.insert(707u32);
    /// hpos.insert(12639u32);
    /// hpos.insert(12638u32);
    /// hpos.insert(818u32);
    /// hpos.insert(2715u32);
    /// let set = HpoSet::new(&ontology, hpos);
    /// assert_eq!(set.information_content().unwrap().gene(), 0.14216587);
    /// ```
    pub fn information_content(&self) -> HpoResult<InformationContent> {
        let n_diseases = self.ontology.omim_diseases().len();
        let n_genes = self.ontology.genes().len();

        let mut ic = InformationContent::default();
        ic.set_gene(n_genes, self.gene_ids().len())?;
        ic.set_omim_disease(n_diseases, self.omim_disease_ids().len())?;
        Ok(ic)
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

    /// Calculates the similarity to another [`HpoSet`]
    ///
    /// # Examples
    ///
    /// ```
    /// use hpo::term::InformationContentKind;
    /// use hpo::{Ontology, HpoSet};
    /// use hpo::term::HpoGroup;
    /// use hpo::similarity::{Builtins, StandardCombiner};
    ///
    /// let ontology = Ontology::from_binary("tests/example.hpo").unwrap();
    ///
    /// // create one set
    /// let mut hpos = HpoGroup::new();
    /// hpos.insert(707u32);
    /// hpos.insert(12639u32);
    /// hpos.insert(12638u32);
    /// hpos.insert(818u32);
    /// hpos.insert(2715u32);
    /// let set = HpoSet::new(&ontology, hpos);
    ///
    /// // create another set
    /// let mut hpos_2 = HpoGroup::new();
    /// hpos_2.insert(100547u32);
    /// hpos_2.insert(12638u32);
    /// hpos_2.insert(864u32);
    /// hpos_2.insert(25454u32);
    /// let set_2 = HpoSet::new(&ontology, hpos_2);
    ///
    /// let similarity = set.similarity(
    ///     &set_2,
    ///     Builtins::GraphIc(InformationContentKind::Omim),
    ///     StandardCombiner::default()
    /// );
    ///
    /// assert_eq!(similarity, 0.8177036);
    /// ```
    pub fn similarity<S: Similarity, C: SimilarityCombiner>(
        &self,
        other: &HpoSet,
        similarity: S,
        combiner: C,
    ) -> f32 {
        let group_sim = GroupSimilarity::new(combiner, similarity);
        group_sim.calculate(self, other)
    }

    /// Returns an Iterator of [`HpoTerm`]
    ///
    /// # Examples
    ///
    /// ```
    /// use hpo::term::InformationContentKind;
    /// use hpo::{Ontology, HpoSet};
    /// use hpo::term::HpoGroup;
    /// use hpo::similarity::{Builtins, StandardCombiner};
    ///
    /// let ontology = Ontology::from_binary("tests/example.hpo").unwrap();
    ///
    /// // create one set
    /// let mut hpos = HpoGroup::new();
    /// hpos.insert(707u32);
    /// hpos.insert(12639u32);
    ///
    /// let set = HpoSet::new(&ontology, hpos);
    ///
    /// let mut iter = set.iter();
    /// assert!(iter.next().is_some());
    /// assert!(iter.next().is_some());
    /// assert!(iter.next().is_none());
    /// ```
    pub fn iter(&self) -> Iter<'_> {
        self.into_iter()
    }

    /// Returns true if the set contains a term with the [`HpoTermId`]
    ///
    /// # Examples
    ///
    /// ```
    /// use hpo::term::InformationContentKind;
    /// use hpo::{Ontology, HpoSet};
    /// use hpo::term::HpoGroup;
    /// use hpo::similarity::{Builtins, StandardCombiner};
    ///
    /// let ontology = Ontology::from_binary("tests/example.hpo").unwrap();
    ///
    /// // create one set
    /// let mut hpos = HpoGroup::new();
    /// hpos.insert(707u32);
    /// hpos.insert(12639u32);
    ///
    /// let set = HpoSet::new(&ontology, hpos);
    ///
    /// assert!(set.contains(&707u32.into()));
    /// assert!(!set.contains(&66666u32.into()));
    /// ```
    pub fn contains(&self, id: &HpoTermId) -> bool {
        self.group.contains(id)
    }
}

impl<'a> IntoIterator for &'a HpoSet<'a> {
    type Item = HpoTerm<'a>;
    type IntoIter = Iter<'a>;
    fn into_iter(self) -> Self::IntoIter {
        Iter::new(self.group.iter(), self.ontology)
    }
}

#[cfg(test)]
mod test {
    use crate::similarity::{Builtins, StandardCombiner};
    use crate::term::internal::HpoTermInternal;
    use crate::term::HpoGroup;
    use crate::term::InformationContentKind;
    use crate::{HpoSet, Ontology};

    #[test]
    fn test() {
        let ontology = Ontology::from_binary("tests/example.hpo").unwrap();

        // create one set
        let mut hpos = HpoGroup::new();
        hpos.insert(707u32);
        hpos.insert(12639u32);
        hpos.insert(12638u32);
        hpos.insert(818u32);
        hpos.insert(2715u32);
        let set = HpoSet::new(&ontology, hpos);
        assert_eq!(set.len(), 5);

        // create another set
        let mut hpos_2 = HpoGroup::new();
        hpos_2.insert(100_547u32);
        hpos_2.insert(12638u32);
        hpos_2.insert(864u32);
        hpos_2.insert(25454u32);
        let set_2 = HpoSet::new(&ontology, hpos_2);
        assert_eq!(set_2.len(), 4);

        let similarity = set.similarity(
            &set_2,
            Builtins::GraphIc(InformationContentKind::Omim),
            StandardCombiner::default(),
        );

        assert!((similarity - 0.817_703_6).abs() < f32::EPSILON);
    }

    #[test]
    fn test_obsolete() {
        let mut ontology = Ontology::from_binary("tests/example.hpo").unwrap();

        let mut obsolete_term = HpoTermInternal::new("Obsolete: Foo".to_string(), 666u32.into());
        *obsolete_term.obsolete_mut() = true;

        ontology.add_term(obsolete_term);

        let mut hpos = HpoGroup::new();
        hpos.insert(707u32);
        hpos.insert(12639u32);
        hpos.insert(12638u32);
        hpos.insert(666u32);
        hpos.insert(2715u32);
        let mut set = HpoSet::new(&ontology, hpos);
        assert_eq!(set.len(), 5);

        let set2 = set.without_obsolete();
        assert_eq!(set.len(), 5);
        assert_eq!(set2.len(), 4);

        set.remove_obsolete();
        assert_eq!(set.len(), 4);
    }

    #[test]
    fn test_with_replaced_obsolete() {
        let mut ontology = Ontology::from_binary("tests/example.hpo").unwrap();

        let mut obsolete_term = HpoTermInternal::new("Obsolete: Foo".to_string(), 666u32.into());
        *obsolete_term.replacement_mut() = Some(25454u32.into());
        ontology.add_term(obsolete_term.clone());

        let mut hpos = HpoGroup::new();
        hpos.insert(707u32);
        hpos.insert(12639u32);
        hpos.insert(12638u32);
        hpos.insert(666u32);
        hpos.insert(2715u32);
        let set = HpoSet::new(&ontology, hpos);
        assert_eq!(set.len(), 5);

        let set2 = set.with_replaced_obsolete();
        assert_eq!(set.len(), 5);
        assert_eq!(set2.len(), 5);

        assert!(set.contains(&666u32.into()));
        assert!(!set2.contains(&666u32.into()));

        assert!(!set.contains(&25454u32.into()));
        assert!(set2.contains(&25454u32.into()));
    }

    #[test]
    fn test_replace_obsolete() {
        let mut ontology = Ontology::from_binary("tests/example.hpo").unwrap();

        let mut obsolete_term = HpoTermInternal::new("Obsolete: Foo".to_string(), 666u32.into());
        *obsolete_term.replacement_mut() = Some(25454u32.into());
        ontology.add_term(obsolete_term.clone());

        let mut hpos = HpoGroup::new();
        hpos.insert(707u32);
        hpos.insert(12639u32);
        hpos.insert(12638u32);
        hpos.insert(666u32);
        hpos.insert(2715u32);
        let mut set = HpoSet::new(&ontology, hpos);
        assert_eq!(set.len(), 5);

        set.replace_obsolete();
        assert_eq!(set.len(), 5);
        assert!(!set.contains(&666u32.into()));
        assert!(set.contains(&25454u32.into()));
    }
}
