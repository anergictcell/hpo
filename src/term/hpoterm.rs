use crate::Similarity;
use crate::annotations::GeneIterator;
use crate::annotations::Genes;
use crate::annotations::OmimDiseaseIterator;
use crate::annotations::OmimDiseases;
use crate::term::HpoTermInternal;
use crate::term::HpoTermIterator;
use crate::HpoParents;
use crate::HpoTermId;
use crate::Ontology;

use crate::HpoError;

use crate::OntologyResult;

use super::HpoGroup;
use super::InformationContent;

#[derive(Debug)]
pub struct HpoTerm<'a> {
    id: &'a HpoTermId,
    name: &'a str,
    parents: &'a HpoParents,
    all_parents: &'a HpoParents,
    genes: &'a Genes,
    omim_diseases: &'a OmimDiseases,
    information_content: &'a InformationContent,
    ontology: &'a Ontology,
}

impl<'a> HpoTerm<'a> {
    pub fn try_new(ontology: &'a Ontology, term: &HpoTermId) -> OntologyResult<HpoTerm<'a>> {
        let term = ontology.get(term).ok_or(HpoError::DoesNotExist)?;
        Ok(HpoTerm {
            id: term.id(),
            name: term.name(),
            parents: term.parents(),
            all_parents: term.all_parents(),
            genes: term.genes(),
            omim_diseases: term.omim_diseases(),
            information_content: term.information_content(),
            ontology,
        })
    }

    pub fn new(ontology: &'a Ontology, term: &'a HpoTermInternal) -> HpoTerm<'a> {
        HpoTerm {
            id: term.id(),
            name: term.name(),
            parents: term.parents(),
            all_parents: term.all_parents(),
            genes: term.genes(),
            omim_diseases: term.omim_diseases(),
            information_content: term.information_content(),
            ontology,
        }
    }

    pub fn id(&self) -> &HpoTermId {
        self.id
    }

    pub fn parents(&self) -> HpoTermIterator<'a> {
        HpoTermIterator::new(self.parents, self.ontology)
    }

    pub fn parent_ids(&self) -> &HpoParents {
        self.parents
    }

    pub fn all_parent_ids(&self) -> &HpoParents {
        self.all_parents
    }

    pub fn all_parents(&self) -> HpoTermIterator<'a> {
        HpoTermIterator::new(self.all_parents, self.ontology)
    }

    pub fn common_ancestor_ids(&self, other: &HpoTerm) -> HpoParents {
        let mut res = self.all_parent_ids() & other.all_parent_ids();

        if other.all_parent_ids().contains(self.id()) {
            res.insert(*self.id());
        }

        if self.all_parent_ids().contains(other.id()) {
            res.insert(*other.id());    
        }

        res
    }

    pub fn union_ancestor_ids(&self, other: &HpoTerm) -> HpoParents {
        self.all_parent_ids() | other.all_parent_ids()
    }

    pub fn common_ancestors(&self, other: &HpoTerm) -> HpoTermOverlap<'a> {
        let group = self.common_ancestor_ids(other);
        HpoTermOverlap::new(group, self.ontology)
    }

    pub fn union_ancestors(&self, other: &HpoTerm) -> HpoTermOverlap<'a> {
        let group = self.union_ancestor_ids(other);
        HpoTermOverlap::new(group, self.ontology)
    }

    pub fn genes(&self) -> GeneIterator<'a> {
        GeneIterator::new(self.genes, self.ontology)
    }

    pub fn omim_diseases(&self) -> OmimDiseaseIterator<'a> {
        OmimDiseaseIterator::new(self.omim_diseases, self.ontology)
    }

    pub fn information_content(&self) -> &InformationContent {
        self.information_content
    }

    pub fn similarity_score(&self, other: &HpoTerm, similarity: &impl Similarity) -> f32 {
        similarity.calculate(self, other)
    }
}

pub struct HpoTermOverlap<'a> {
    overlap: HpoGroup,
    ontology: &'a Ontology,
}

impl<'a> HpoTermOverlap<'a> {
    fn new(overlap: HpoGroup, ontology: &'a Ontology) -> Self {
        Self {
            overlap,
            ontology
        }
    }
}

impl<'a> Iterator for HpoTermOverlap<'a> {
    type Item = HpoTerm<'a>;
    fn next(&mut self) -> Option<Self::Item> {
        match self.overlap.pop() {
            Some(x) => {
                let term = self.ontology.get_unchecked(&x);
                Some(HpoTerm::new(
                    self.ontology,
                    term
                ))
            },
            None => None,
        }
    }
}
