use crate::annotations::GeneIterator;
use crate::annotations::Genes;
use crate::HpoTermId;
use crate::HpoParents;
use crate::Ontology;
use crate::term::HpoTermInternal;
use crate::term::HpoTermIterator;

use crate::HpoError;

use crate::OntologyResult;


#[derive(Debug)]
pub struct HpoTerm<'a> {
    id: &'a HpoTermId,
    name: &'a str,
    parents: &'a HpoParents,
    all_parents: &'a HpoParents,
    genes: &'a Genes,
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

    pub fn overlap(&self, other: &HpoTerm) -> HpoParents {
        self.all_parents & other.all_parents
    }

    pub fn genes(&self) -> GeneIterator<'a> {
        GeneIterator::new(self.genes, self.ontology)
    }
}
