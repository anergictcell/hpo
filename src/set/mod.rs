//! An `HpoSet` can represent e.g. the clinical information of a patient or the symptoms of a disease
use crate::HpoTerm;
use crate::annotations::Genes;
use crate::annotations::OmimDiseases;
use crate::term::HpoGroup;
use crate::term::HpoTerms;
use crate::Ontology;
use crate::term::InformationContent;

/// A set of unique HPO terms
///
/// As in a set, each term can only appear once
/// though that is not yet guaranteed in the implementation (TODO)
pub struct HpoSet<'a> {
    ontology: &'a Ontology,
    group: HpoGroup
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
        let group = self.group
            .iter()
            .filter(|term1_id| {
                self.group.iter().all(|term2_id| {
                    self.ontology.get_unchecked(term2_id).all_parents().contains(term1_id)
                })
            })
            .collect();
            HpoSet { ontology: self.ontology, group }
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

    /// Returns all [`Genes`] that are associated to the set
    ///
    /// # Panics
    ///
    /// When an `HpoTermId` of the set is not part of the Ontology
    fn gene_ids(&self) -> Genes {
        self.group.into_iter()
        .map(|term_id| self.ontology.get_unchecked(term_id).genes())
        .fold(Genes::default(), |acc, element| {&acc | element})
    }

    /// Returns all [`OmimDiseases`] that are associated to the set
    ///
    /// # Panics
    ///
    /// When an `HpoTermId` of the set is not part of the Ontology
    fn omim_disease_ids(&self) -> OmimDiseases {
        self.group.into_iter()
        .map(|term_id| self.ontology.get_unchecked(term_id).omim_diseases())
        .fold(OmimDiseases::default(), |acc, element| {&acc | element})
    }

    /// Calculates and returns the aggregated [`InformationContent`] of the set
    ///
    /// The `InformationContent` is not cached internally, so this operation
    /// is not cheap
    ///
    /// # Panics
    ///
    /// - When an `HpoTermId` of the set is not part of the Ontology
    /// - When the ontology or set have more than u16::MAX genes or diseases
    pub fn information_content(&self) -> InformationContent {
        let n_diseases = self.ontology.omim_diseases().len();
        let n_genes = self.ontology.genes().len();

        let mut ic = InformationContent::default();
        ic.set_gene(n_genes, self.gene_ids().len()).expect("Too many genes caused overflow");
        ic.set_omim_disease(n_diseases, self.omim_disease_ids().len()).expect("Too many genes caused overflow");
        ic
    }

    pub fn common_ancestor_ids(&self) -> HpoGroup {
        unimplemented!()
    }

    pub fn combinations(&self) -> BothWayCombinations {
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
}

impl<'a> IntoIterator for &'a HpoSet<'a> {
    type Item = HpoTerm<'a>;
    type IntoIter = HpoTerms<'a>;
    fn into_iter(self) -> Self::IntoIter {
        HpoTerms::new(&self.group, self.ontology)
    }
}

// impl<'a> From<HpoTerms<'a>> for HpoSet<'a> {
//     fn from(terms: HpoTerms<'a>) -> HpoSet<'a> {
//         let ontology = terms.ontology();
//         let group = HpoGroup::default();
//         for term in &terms {
//             group.insert(term.id());
//         }
//         // = &terms.map(|t| *t.id()).collect();
//         Self {
//             ontology,
//             group,
//         }
//     }
// }

pub struct BothWayCombinations<'a> {
    group: HpoSet<'a>,
    index_a: usize,
    index_b: usize
}

impl<'a> Iterator for BothWayCombinations<'a> {
    type Item = (HpoTerm<'a>, HpoTerm<'a>);
    fn next(&mut self) -> Option<Self::Item> {
        if self.index_a > self.group.len() {
            return None
        }

        let b = if let Some(b) = self.group.get(self.index_b) {
            b
        } else {
            self.index_a += 1;
            self.index_b = 0;
            return self.next()
        };

        let a = self.group.get(self.index_a)?;

        self.index_b += 1;

        Some((a, b))
    }
}

pub struct OneWayCombinations<'a> {
    group: HpoSet<'a>,
    index_a: usize,
    index_b: usize
}

impl<'a> Iterator for OneWayCombinations<'a> {
    type Item = (HpoTerm<'a>, HpoTerm<'a>);
    fn next(&mut self) -> Option<Self::Item> {
        if self.index_a > self.group.len() {
            return None
        }

        let b = if let Some(b) = self.group.get(self.index_b) {
            b
        } else {
            self.index_a += 1;
            self.index_b = self.index_a + 1;
            return self.next()
        };

        let a = self.group.get(self.index_a)?;

        self.index_b += 1;

        Some((a, b))
    }
}