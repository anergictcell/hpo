use crate::annotations::GeneIterator;
use crate::annotations::Genes;
use crate::annotations::OmimDiseaseIterator;
use crate::annotations::OmimDiseases;
use crate::similarity::Similarity;
use crate::term::internal::HpoTermInternal;
use crate::term::HpoParents;
use crate::term::HpoTerms;
use crate::HpoTermId;
use crate::Ontology;

use crate::HpoError;

use crate::HpoResult;

use super::group::GroupCombine;
use super::HpoChildren;
use super::InformationContent;

/// The `HpoTerm` represents a single term from the HP Ontology
///
/// The term holds all required information and relationship data.
/// It provides functionality for path traversals and similarity calculations.
#[derive(Debug, Clone, Copy)]
pub struct HpoTerm<'a> {
    id: &'a HpoTermId,
    name: &'a str,
    parents: &'a HpoParents,
    all_parents: &'a HpoParents,
    children: &'a HpoChildren,
    genes: &'a Genes,
    omim_diseases: &'a OmimDiseases,
    information_content: &'a InformationContent,
    ontology: &'a Ontology,
}

impl<'a> HpoTerm<'a> {
    /// Constructs a new [`HpoTerm`]
    ///
    /// # Errors
    ///
    /// If the given [`HpoTermId`] does not match an existing term
    /// it returns an Error
    pub fn try_new(ontology: &'a Ontology, term: &HpoTermId) -> HpoResult<HpoTerm<'a>> {
        let term = ontology.get(term).ok_or(HpoError::DoesNotExist)?;
        Ok(HpoTerm::new(ontology, term))
    }

    /// Constructs a new [`HpoTerm`] from an `HpoTermInternal`
    pub(crate) fn new(ontology: &'a Ontology, term: &'a HpoTermInternal) -> HpoTerm<'a> {
        HpoTerm {
            id: term.id(),
            name: term.name(),
            parents: term.parents(),
            all_parents: term.all_parents(),
            children: term.children(),
            genes: term.genes(),
            omim_diseases: term.omim_diseases(),
            information_content: term.information_content(),
            ontology,
        }
    }

    /// Returns the [`HpoTermId`] of the term
    ///
    /// e.g.: `HP:0012345`
    pub fn id(&self) -> &HpoTermId {
        self.id
    }

    /// Returns the name of the term
    ///
    /// e.g.: `Abnormality of the nervous system`
    pub fn name(&self) -> &str {
        self.name
    }

    /// Returns an iterator of the direct patients of the term
    pub fn parents(&self) -> HpoTerms<'a> {
        HpoTerms::new(self.parents, self.ontology)
    }

    /// Returns an iterator of the direct children of the term
    pub fn children(&self) -> HpoTerms<'a> {
        HpoTerms::new(self.children, self.ontology)
    }

    /// Returns the [`HpoTermId`]s of the direct parents
    pub fn parent_ids(&self) -> &HpoParents {
        self.parents
    }

    /// Returns the [`HpoTermId`]s of al; direct and indirect parents
    pub fn all_parent_ids(&self) -> &HpoParents {
        self.all_parents
    }

    /// Returns an iterator of the direct and indrect patients of the term
    pub fn all_parents(&self) -> HpoTerms<'a> {
        HpoTerms::new(self.all_parents, self.ontology)
    }

    /// Returns the [`HpoTermId`]s that are parents of both `self` **and** `other`
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

    /// Returns the [`HpoTermId`]s that are parents of either `self` **or** `other`
    pub fn union_ancestor_ids(&self, other: &HpoTerm) -> HpoParents {
        self.all_parent_ids() | other.all_parent_ids()
    }

    /// Returns an iterator of [`HpoTerm`]s that are parents of both `self` **and** `other`
    pub fn common_ancestors(&self, other: &HpoTerm) -> GroupCombine {
        GroupCombine::new(self.common_ancestor_ids(other), self.ontology)
    }

    /// Returns an iterator of [`HpoTerm`]s that are parents of either `self` **or** `other`
    pub fn union_ancestors(&self, other: &HpoTerm) -> GroupCombine {
        GroupCombine::new(self.union_ancestor_ids(other), self.ontology)
    }

    /// Returns an iterator of all associated [`crate::annotations::Gene`]s
    pub fn genes(&self) -> GeneIterator<'a> {
        GeneIterator::new(self.genes, self.ontology)
    }

    /// Returns an iterator of all associated [`crate::annotations::OmimDisease`]s
    pub fn omim_diseases(&self) -> OmimDiseaseIterator<'a> {
        OmimDiseaseIterator::new(self.omim_diseases, self.ontology)
    }

    /// Returns the [`InformationContent`] of the term
    pub fn information_content(&self) -> &InformationContent {
        self.information_content
    }

    /// Calculates the similarity of `self` and `other` using the provided `Similarity` algorithm
    pub fn similarity_score(&self, other: &HpoTerm, similarity: &impl Similarity) -> f32 {
        similarity.calculate(self, other)
    }

    /// Returns the distance (steps) from `self` to `other`, if `other` is a parent of `self`
    pub fn distance_to_ancestor(&self, other: &HpoTerm) -> Option<usize> {
        if self.id() == other.id() {
            return Some(0);
        }
        if self.parent_ids().contains(other.id()) {
            return Some(1);
        }
        if !self.all_parent_ids().contains(other.id()) {
            return None;
        }
        self.parents()
            .filter_map(|p| p.distance_to_ancestor(other))
            .min()
            .map(|c| c + 1)
    }

    /// Returns `true` if `self` is a child (direct or indirect) of `other`
    pub fn child_of(&self, other: &HpoTerm) -> bool {
        self.all_parent_ids().contains(other.id())
    }

    /// Returns `true` if `self` is a parent (direct or indirect) of `other`
    pub fn parent_of(&self, other: &HpoTerm) -> bool {
        other.child_of(self)
    }

    /// Returns the shortest path to traverse from `self` ot `other`, if `other` is a parent of `self`
    pub fn path_to_ancestor(&self, other: &HpoTerm) -> Option<Vec<HpoTermId>> {
        if self.id() == other.id() {
            return Some(vec![]);
        }
        if self.parent_ids().contains(other.id()) {
            return Some(vec![*other.id()]);
        }
        if !self.all_parent_ids().contains(other.id()) {
            return None;
        }
        self.parents()
            .filter_map(|p| match p.path_to_ancestor(other) {
                Some(mut x) => {
                    x.insert(0, *p.id);
                    Some(x)
                }
                None => None,
            })
            .min_by_key(Vec::len)
    }

    /// Returns the distance (steps) from `self` to `other`
    pub fn distance_to_term(&self, other: &HpoTerm) -> Option<usize> {
        self.common_ancestors(other)
            .filter_map(|parent| {
                Some(self.distance_to_ancestor(&parent)? + other.distance_to_ancestor(&parent)?)
            })
            .min()
    }
}
