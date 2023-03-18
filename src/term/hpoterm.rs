use crate::annotations::GeneIterator;
use crate::annotations::Genes;
use crate::annotations::OmimDiseaseIterator;
use crate::annotations::OmimDiseases;
use crate::similarity::Similarity;
use crate::term::internal::HpoTermInternal;
use crate::term::HpoGroup;
use crate::term::Iter;
use crate::HpoTermId;
use crate::Ontology;

use crate::HpoError;

use crate::HpoResult;

use super::group::Combined;

use super::InformationContent;

/// The `HpoTerm` represents a single term from the HP Ontology
///
/// The term holds all required information and relationship data.
/// It provides functionality for path traversals and similarity calculations.
#[derive(Debug, Clone, Copy)]
pub struct HpoTerm<'a> {
    id: &'a HpoTermId,
    name: &'a str,
    parents: &'a HpoGroup,
    all_parents: &'a HpoGroup,
    children: &'a HpoGroup,
    genes: &'a Genes,
    omim_diseases: &'a OmimDiseases,
    information_content: &'a InformationContent,
    obsolete: bool,
    replaced_by: Option<HpoTermId>,
    ontology: &'a Ontology,
}

impl<'a> HpoTerm<'a> {
    /// Constructs a new [`HpoTerm`]
    ///
    /// # Errors
    ///
    /// If the given [`HpoTermId`] does not match an existing term
    /// it returns an Error
    ///
    /// # Examples
    ///
    /// ```
    /// use hpo::{HpoTerm, Ontology};
    ///
    /// let ontology = Ontology::from_binary("tests/example.hpo").unwrap();
    ///
    /// let term = HpoTerm::try_new(&ontology, 118u32.into());
    /// assert!(term.is_ok());
    ///
    /// let non_existing_term = HpoTerm::try_new(&ontology, 666666666u32.into());
    /// assert!(non_existing_term.is_err());
    /// ```
    pub fn try_new(ontology: &'a Ontology, term: HpoTermId) -> HpoResult<HpoTerm<'a>> {
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
            obsolete: term.obsolete(),
            replaced_by: term.replacement(),
            ontology,
        }
    }

    /// Returns the [`HpoTermId`] of the term
    ///
    /// e.g.: `HP:0012345`
    ///
    /// # Examples
    ///
    /// ```
    /// use hpo::{HpoTerm, Ontology};
    ///
    /// let ontology = Ontology::from_binary("tests/example.hpo").unwrap();
    ///
    /// let term = HpoTerm::try_new(&ontology, 118u32.into()).unwrap();
    /// assert_eq!(term.id(), "HP:0000118");
    /// ```
    pub fn id(&self) -> HpoTermId {
        *self.id
    }

    /// Returns the name of the term
    ///
    /// e.g.: `Abnormality of the nervous system`
    ///
    /// # Examples
    ///
    /// ```
    /// use hpo::{HpoTerm, Ontology};
    ///
    /// let ontology = Ontology::from_binary("tests/example.hpo").unwrap();
    ///
    /// let term = HpoTerm::try_new(&ontology, 118u32.into()).unwrap();
    /// assert_eq!(term.name(), "Phenotypic abnormality");
    /// ```
    pub fn name(&self) -> &str {
        self.name
    }

    /// Returns the [`HpoTermId`]s of the direct parents
    ///
    /// # Examples
    ///
    /// ```
    /// use hpo::{HpoTerm, Ontology};
    ///
    /// let ontology = Ontology::from_binary("tests/example.hpo").unwrap();
    ///
    /// let term = HpoTerm::try_new(&ontology, 11017u32.into()).unwrap();
    /// assert_eq!(term.parent_ids().len(), 1);
    /// ```
    pub fn parent_ids(&self) -> &HpoGroup {
        self.parents
    }

    /// Returns an iterator of the direct patients of the term
    ///
    /// # Examples
    ///
    /// ```
    /// use hpo::{HpoTerm, Ontology};
    ///
    /// let ontology = Ontology::from_binary("tests/example.hpo").unwrap();
    ///
    /// let term = HpoTerm::try_new(&ontology, 11017u32.into()).unwrap();
    /// assert_eq!(term.parents().count(), 1);
    /// for parent in term.parents() {
    ///    println!("{}", parent.name());
    /// }
    /// ```
    pub fn parents(&self) -> Iter<'a> {
        Iter::new(self.parents.iter(), self.ontology)
    }

    /// Returns the [`HpoTermId`]s of al; direct and indirect parents
    ///
    /// # Examples
    ///
    /// ```
    /// use hpo::{HpoTerm, Ontology};
    ///
    /// let ontology = Ontology::from_binary("tests/example.hpo").unwrap();
    ///
    /// let term = HpoTerm::try_new(&ontology, 11017u32.into()).unwrap();
    /// assert_eq!(term.all_parent_ids().len(), 3);
    /// ```
    pub fn all_parent_ids(&self) -> &HpoGroup {
        self.all_parents
    }

    /// Returns an iterator of the direct and indrect patients of the term
    ///
    /// # Examples
    ///
    /// ```
    /// use hpo::{HpoTerm, Ontology};
    ///
    /// let ontology = Ontology::from_binary("tests/example.hpo").unwrap();
    ///
    /// let term = HpoTerm::try_new(&ontology, 11017u32.into()).unwrap();
    /// assert_eq!(term.all_parents().count(), 3);
    /// for parent in term.all_parents() {
    ///    println!("{}", parent.name());
    /// }
    /// ```
    pub fn all_parents(&self) -> Iter<'a> {
        Iter::new(self.all_parents.iter(), self.ontology)
    }

    /// Returns the [`HpoTermId`]s of the direct children
    ///
    /// # Examples
    ///
    /// ```
    /// use hpo::{HpoTerm, Ontology};
    ///
    /// let ontology = Ontology::from_binary("tests/example.hpo").unwrap();
    ///
    /// let term = HpoTerm::try_new(&ontology, 1939u32.into()).unwrap();
    /// assert_eq!(term.children_ids().len(), 2);
    /// ```
    pub fn children_ids(&self) -> &HpoGroup {
        self.children
    }

    /// Returns an iterator of the direct children of the term
    ///
    /// # Examples
    ///
    /// ```
    /// use hpo::{HpoTerm, Ontology};
    ///
    /// let ontology = Ontology::from_binary("tests/example.hpo").unwrap();
    ///
    /// let term = HpoTerm::try_new(&ontology, 1939u32.into()).unwrap();
    /// for child in term.children() {
    ///    println!("{}", child.name());
    /// }
    /// ```
    pub fn children(&self) -> Iter<'a> {
        Iter::new(self.children.iter(), self.ontology)
    }

    /// Returns the [`HpoTermId`]s that are parents of both `self` **and** `other`
    ///
    /// # Note:
    ///
    /// This method includes `self` and `other` into their corresponding
    /// parent-id list so that if one term is a parent term of the other,
    /// it is included. It also means that `self` is included if `self == other`.
    ///
    /// This might seem counter-intuitive at first, but in most cases this
    /// is what is actually needed
    ///
    /// # Examples
    ///
    /// ## Note:
    /// Check out the [`crate::Ontology`] documentation to see the example ontology
    ///
    /// ## Test-case 1: Terms do **not** have a  parent - child relationship
    /// ```
    /// use hpo::{HpoTerm, Ontology};
    ///
    /// let ontology = Ontology::from_binary("tests/example.hpo").unwrap();
    ///
    /// let term1 = HpoTerm::try_new(&ontology, 11017u32.into()).unwrap();
    /// let term2 = HpoTerm::try_new(&ontology, 25454u32.into()).unwrap();
    ///
    /// assert!(!term1.parent_of(&term2));
    /// assert!(!term2.parent_of(&term1));
    ///
    /// let common_ancestors = term1.all_common_ancestor_ids(&term2);
    /// assert_eq!(common_ancestors.len(), 3);
    ///
    /// assert!(!common_ancestors.contains(&11017u32.into()));
    /// assert!(!common_ancestors.contains(&25454u32.into()));
    /// ```
    ///
    /// ## Test-case 2: Terms **do have** a  parent - child relationship
    /// ```
    /// use hpo::{HpoTerm, Ontology};
    ///
    /// let ontology = Ontology::from_binary("tests/example.hpo").unwrap();
    ///
    /// let term1 = HpoTerm::try_new(&ontology, 11017u32.into()).unwrap();
    /// let term2 = HpoTerm::try_new(&ontology, 1939u32.into()).unwrap();
    ///
    /// assert!(!term1.parent_of(&term2));
    /// assert!(term2.parent_of(&term1));
    ///
    /// let common_ancestors = term1.all_common_ancestor_ids(&term2);
    /// assert_eq!(common_ancestors.len(), 3);
    ///
    /// assert!(!common_ancestors.contains(&11017u32.into()));
    /// assert!(common_ancestors.contains(&1939u32.into()));
    /// ```
    pub fn all_common_ancestor_ids(&self, other: &HpoTerm) -> HpoGroup {
        (self.all_parent_ids() + self.id()) & (other.all_parent_ids() + other.id())
    }

    /// Returns the [`HpoTermId`]s that are parents of both `self` **and** `other`
    ///
    /// # Note:
    ///
    /// This method does not include `self` and `other`. Depending on your
    /// use case, you might prefer [`HpoTerm::all_common_ancestor_ids`].
    ///
    /// # Examples
    ///
    /// ## Note:
    /// Check out the [`crate::Ontology`] documentation to see the example ontology
    ///
    /// ## Test-case 1: Terms do **not** have a  parent - child relationship
    /// ```
    /// use hpo::{HpoTerm, Ontology};
    ///
    /// let ontology = Ontology::from_binary("tests/example.hpo").unwrap();
    ///
    /// let term1 = HpoTerm::try_new(&ontology, 11017u32.into()).unwrap();
    /// let term2 = HpoTerm::try_new(&ontology, 25454u32.into()).unwrap();
    ///
    /// assert!(!term1.parent_of(&term2));
    /// assert!(!term2.parent_of(&term1));
    ///
    /// let common_ancestors = term1.common_ancestor_ids(&term2);
    /// assert_eq!(common_ancestors.len(), 3);
    ///
    /// assert!(!common_ancestors.contains(&11017u32.into()));
    /// assert!(!common_ancestors.contains(&25454u32.into()));
    /// ```
    ///
    /// ## Test-case 2: Terms **do have** a  parent - child relationship
    /// ```
    /// use hpo::{HpoTerm, Ontology};
    ///
    /// let ontology = Ontology::from_binary("tests/example.hpo").unwrap();
    ///
    /// let term1 = HpoTerm::try_new(&ontology, 11017u32.into()).unwrap();
    /// let term2 = HpoTerm::try_new(&ontology, 1939u32.into()).unwrap();
    ///
    /// assert!(!term1.parent_of(&term2));
    /// assert!(term2.parent_of(&term1));
    ///
    /// let common_ancestors = term1.common_ancestor_ids(&term2);
    /// assert_eq!(common_ancestors.len(), 2);
    ///
    /// assert!(!common_ancestors.contains(&11017u32.into()));
    /// assert!(!common_ancestors.contains(&1939u32.into()));
    /// ```
    pub fn common_ancestor_ids(&self, other: &HpoTerm) -> HpoGroup {
        self.all_parent_ids() & other.all_parent_ids()
    }

    /// Returns the [`HpoTermId`]s that are parents of either `self` **or** `other`
    ///
    /// # Note:
    ///
    /// This method includes `self` and `other` into their corresponding
    /// parent-id list so that both are included themselves as well.
    ///
    /// This might seem counter-intuitive at first, but in many cases this
    /// is what is actually needed
    ///
    /// # Examples
    ///
    /// ## Note:
    /// Check out the [`crate::Ontology`] documentation to see the example ontology
    ///
    /// ## Test-case 1: Terms do **not** have a  parent - child relationship
    /// ```
    /// use hpo::{HpoTerm, Ontology};
    ///
    /// let ontology = Ontology::from_binary("tests/example.hpo").unwrap();
    ///
    /// let term1 = HpoTerm::try_new(&ontology, 11017u32.into()).unwrap();
    /// let term2 = HpoTerm::try_new(&ontology, 12639u32.into()).unwrap();
    ///
    /// assert!(!term1.parent_of(&term2));
    /// assert!(!term2.parent_of(&term1));
    ///
    /// let union_ancestors = term1.all_union_ancestor_ids(&term2);
    /// assert_eq!(union_ancestors.len(), 4);
    ///
    /// assert!(!union_ancestors.contains(&11017u32.into()));
    /// assert!(!union_ancestors.contains(&12639u32.into()));
    /// ```
    ///
    /// ## Test-case 2: Terms **do have** a  parent - child relationship
    /// ```
    /// use hpo::{HpoTerm, Ontology};
    ///
    /// let ontology = Ontology::from_binary("tests/example.hpo").unwrap();
    ///
    /// let term1 = HpoTerm::try_new(&ontology, 11017u32.into()).unwrap();
    /// let term2 = HpoTerm::try_new(&ontology, 1939u32.into()).unwrap();
    ///
    /// assert!(!term1.parent_of(&term2));
    /// assert!(term2.parent_of(&term1));
    ///
    /// let union_ancestors = term1.all_union_ancestor_ids(&term2);
    /// assert_eq!(union_ancestors.len(), 3);
    ///
    /// assert!(!union_ancestors.contains(&11017u32.into()));
    /// assert!(union_ancestors.contains(&1939u32.into()));
    /// ```
    pub fn all_union_ancestor_ids(&self, other: &HpoTerm) -> HpoGroup {
        self.all_parent_ids() | other.all_parent_ids()
    }

    /// Returns the [`HpoTermId`]s that are parents of either `self` **or** `other`
    ///
    /// # Examples
    ///
    /// ## Note:
    /// Check out the [`crate::Ontology`] documentation to see the example ontology
    ///
    /// ## Test-case 1: Terms do **not** have a  parent - child relationship
    /// ```
    /// use hpo::{HpoTerm, Ontology};
    ///
    /// let ontology = Ontology::from_binary("tests/example.hpo").unwrap();
    ///
    /// let term1 = HpoTerm::try_new(&ontology, 11017u32.into()).unwrap();
    /// let term2 = HpoTerm::try_new(&ontology, 12639u32.into()).unwrap();
    ///
    /// assert!(!term1.parent_of(&term2));
    /// assert!(!term2.parent_of(&term1));
    ///
    /// let union_ancestors = term1.all_union_ancestor_ids(&term2);
    /// assert_eq!(union_ancestors.len(), 4);
    ///
    /// assert!(!union_ancestors.contains(&11017u32.into()));
    /// assert!(!union_ancestors.contains(&12639u32.into()));
    /// ```
    ///
    /// ## Test-case 2: Terms **do have** a  parent - child relationship
    /// ```
    /// use hpo::{HpoTerm, Ontology};
    ///
    /// let ontology = Ontology::from_binary("tests/example.hpo").unwrap();
    ///
    /// let term1 = HpoTerm::try_new(&ontology, 11017u32.into()).unwrap();
    /// let term2 = HpoTerm::try_new(&ontology, 1939u32.into()).unwrap();
    ///
    /// assert!(!term1.parent_of(&term2));
    /// assert!(term2.parent_of(&term1));
    ///
    /// let union_ancestors = term1.all_union_ancestor_ids(&term2);
    /// assert_eq!(union_ancestors.len(), 3);
    ///
    /// assert!(!union_ancestors.contains(&11017u32.into()));
    /// assert!(union_ancestors.contains(&1939u32.into()));
    /// ```
    pub fn union_ancestor_ids(&self, other: &HpoTerm) -> HpoGroup {
        self.all_parent_ids() | other.all_parent_ids()
    }

    /// Returns an iterator of [`HpoTerm`]s that are parents of both `self` **and** `other`
    ///
    /// # Note:
    ///
    /// This method includes `self` and `other` into their corresponding
    /// parent-id list so that if one term is a parent term of the other,
    /// it is included. It also means that `self` is included if `self == other`.
    ///
    /// This might seem counter-intuitive at first, but in most cases this
    /// is what is actually needed
    ///
    /// # Examples
    ///
    /// ## Note:
    /// Check out the [`crate::Ontology`] documentation to see the example ontology
    ///
    /// ## Test-case 1: Terms do **not** have a  parent - child relationship
    /// ```
    /// use hpo::{HpoTerm, Ontology};
    ///
    /// let ontology = Ontology::from_binary("tests/example.hpo").unwrap();
    ///
    /// let term1 = HpoTerm::try_new(&ontology, 11017u32.into()).unwrap();
    /// let term2 = HpoTerm::try_new(&ontology, 25454u32.into()).unwrap();
    ///
    /// assert!(!term1.parent_of(&term2));
    /// assert!(!term2.parent_of(&term1));
    ///
    /// let common_ancestors = term1.all_common_ancestors(&term2);
    /// assert_eq!(common_ancestors.len(), 3);
    ///
    /// for term in &common_ancestors {
    ///     println!("{}", term.name());
    /// }
    /// ```
    ///
    /// ## Test-case 2: Terms **do have** a  parent - child relationship
    /// ```
    /// use hpo::{HpoTerm, Ontology};
    ///
    /// let ontology = Ontology::from_binary("tests/example.hpo").unwrap();
    ///
    /// let term1 = HpoTerm::try_new(&ontology, 11017u32.into()).unwrap();
    /// let term2 = HpoTerm::try_new(&ontology, 1939u32.into()).unwrap();
    ///
    /// assert!(!term1.parent_of(&term2));
    /// assert!(term2.parent_of(&term1));
    ///
    /// let common_ancestors = term1.all_common_ancestors(&term2);
    /// assert_eq!(common_ancestors.len(), 3);
    ///
    /// for term in &common_ancestors {
    ///     println!("{}", term.name());
    /// }
    /// ```
    pub fn all_common_ancestors(&self, other: &HpoTerm) -> Combined {
        Combined::new(self.all_common_ancestor_ids(other), self.ontology)
    }

    /// Returns an iterator of [`HpoTerm`]s that are parents of both `self` **and** `other`
    ///
    /// # Note:
    ///
    /// This method does not include `self` and `other`. Depending on your
    /// use case, you might prefer [`HpoTerm.all_common_ancestor_ids`].
    ///
    /// ## Test-case 1: Terms do **not** have a  parent - child relationship
    /// ```
    /// use hpo::{HpoTerm, Ontology};
    ///
    /// let ontology = Ontology::from_binary("tests/example.hpo").unwrap();
    ///
    /// let term1 = HpoTerm::try_new(&ontology, 11017u32.into()).unwrap();
    /// let term2 = HpoTerm::try_new(&ontology, 25454u32.into()).unwrap();
    ///
    /// assert!(!term1.parent_of(&term2));
    /// assert!(!term2.parent_of(&term1));
    ///
    /// let common_ancestors = term1.common_ancestors(&term2);
    /// assert_eq!(common_ancestors.len(), 3);
    ///
    /// for term in &common_ancestors {
    ///     println!("{}", term.name());
    /// }
    /// ```
    ///
    /// ## Test-case 2: Terms **do have** a  parent - child relationship
    /// ```
    /// use hpo::{HpoTerm, Ontology};
    ///
    /// let ontology = Ontology::from_binary("tests/example.hpo").unwrap();
    ///
    /// let term1 = HpoTerm::try_new(&ontology, 11017u32.into()).unwrap();
    /// let term2 = HpoTerm::try_new(&ontology, 1939u32.into()).unwrap();
    ///
    /// assert!(!term1.parent_of(&term2));
    /// assert!(term2.parent_of(&term1));
    ///
    /// let common_ancestors = term1.common_ancestors(&term2);
    /// assert_eq!(common_ancestors.len(), 2);
    ///
    /// for term in &common_ancestors {
    ///     println!("{}", term.name());
    /// }
    /// ```
    pub fn common_ancestors(&self, other: &HpoTerm) -> Combined {
        Combined::new(self.common_ancestor_ids(other), self.ontology)
    }

    /// Returns an iterator of [`HpoTerm`]s that are parents of either `self` **or** `other`
    ///
    /// # Note:
    ///
    /// This method includes `self` and `other` into their corresponding
    /// parent-id list so that both are included themselves as well.
    ///
    /// This might seem counter-intuitive at first, but in many cases this
    /// is what is actually needed
    ///
    /// # Examples
    ///
    /// ## Note:
    /// Check out the [`crate::Ontology`] documentation to see the example ontology
    ///
    /// ## Test-case 1: Terms do **not** have a  parent - child relationship
    /// ```
    /// use hpo::{HpoTerm, Ontology};
    ///
    /// let ontology = Ontology::from_binary("tests/example.hpo").unwrap();
    ///
    /// let term1 = HpoTerm::try_new(&ontology, 11017u32.into()).unwrap();
    /// let term2 = HpoTerm::try_new(&ontology, 12639u32.into()).unwrap();
    ///
    /// assert!(!term1.parent_of(&term2));
    /// assert!(!term2.parent_of(&term1));
    ///
    /// let union_ancestors = term1.all_union_ancestors(&term2);
    /// assert_eq!(union_ancestors.len(), 4);
    ///
    /// for term in &union_ancestors {
    ///     println!("{}", term.name());
    /// }
    /// ```
    ///
    /// ## Test-case 2: Terms **do have** a  parent - child relationship
    /// ```
    /// use hpo::{HpoTerm, Ontology};
    ///
    /// let ontology = Ontology::from_binary("tests/example.hpo").unwrap();
    ///
    /// let term1 = HpoTerm::try_new(&ontology, 11017u32.into()).unwrap();
    /// let term2 = HpoTerm::try_new(&ontology, 1939u32.into()).unwrap();
    ///
    /// assert!(!term1.parent_of(&term2));
    /// assert!(term2.parent_of(&term1));
    ///
    /// let union_ancestors = term1.all_union_ancestors(&term2);
    /// assert_eq!(union_ancestors.len(), 3);
    ///
    /// for term in &union_ancestors {
    ///     println!("{}", term.name());
    /// }
    /// ```
    pub fn all_union_ancestors(&self, other: &HpoTerm) -> Combined {
        Combined::new(self.union_ancestor_ids(other), self.ontology)
    }

    /// Returns an iterator of [`HpoTerm`]s that are parents of either `self` **or** `other`
    ///
    /// # Examples
    ///
    /// ## Note:
    /// Check out the [`crate::Ontology`] documentation to see the example ontology
    ///
    /// ## Test-case 1: Terms do **not** have a  parent - child relationship
    /// ```
    /// use hpo::{HpoTerm, Ontology};
    ///
    /// let ontology = Ontology::from_binary("tests/example.hpo").unwrap();
    ///
    /// let term1 = HpoTerm::try_new(&ontology, 11017u32.into()).unwrap();
    /// let term2 = HpoTerm::try_new(&ontology, 12639u32.into()).unwrap();
    ///
    /// assert!(!term1.parent_of(&term2));
    /// assert!(!term2.parent_of(&term1));
    ///
    /// let union_ancestors = term1.union_ancestors(&term2);
    /// assert_eq!(union_ancestors.len(), 4);
    ///
    /// for term in &union_ancestors {
    ///     println!("{}", term.name());
    /// }
    /// ```
    ///
    /// ## Test-case 2: Terms **do have** a  parent - child relationship
    /// ```
    /// use hpo::{HpoTerm, Ontology};
    ///
    /// let ontology = Ontology::from_binary("tests/example.hpo").unwrap();
    ///
    /// let term1 = HpoTerm::try_new(&ontology, 11017u32.into()).unwrap();
    /// let term2 = HpoTerm::try_new(&ontology, 1939u32.into()).unwrap();
    ///
    /// assert!(!term1.parent_of(&term2));
    /// assert!(term2.parent_of(&term1));
    ///
    /// let union_ancestors = term1.union_ancestors(&term2);
    /// assert_eq!(union_ancestors.len(), 3);
    ///
    /// for term in &union_ancestors {
    ///     println!("{}", term.name());
    /// }
    /// ```
    pub fn union_ancestors(&self, other: &HpoTerm) -> Combined {
        Combined::new(self.union_ancestor_ids(other), self.ontology)
    }

    /// Returns an iterator of all associated [`crate::annotations::Gene`]s
    ///
    /// # Examples
    ///
    /// ```
    /// use hpo::{HpoTerm, Ontology};
    ///
    /// let ontology = Ontology::from_binary("tests/example.hpo").unwrap();
    ///
    /// let term = HpoTerm::try_new(&ontology, 11017u32.into()).unwrap();
    /// for gene in term.genes() {
    ///     println!("{}", gene.name());
    /// }
    /// ```
    pub fn genes(&self) -> GeneIterator<'a> {
        GeneIterator::new(self.genes, self.ontology)
    }

    /// Returns the set of all associated [`crate::annotations::Gene`]s
    ///
    /// # Examples
    ///
    /// ```
    /// use hpo::{HpoTerm, Ontology};
    ///
    /// let ontology = Ontology::from_binary("tests/example.hpo").unwrap();
    ///
    /// let term = HpoTerm::try_new(&ontology, 11017u32.into()).unwrap();
    /// assert_eq!(term.gene_ids().len(), 575);
    /// ```
    pub fn gene_ids(&self) -> &Genes {
        self.genes
    }

    /// Returns an iterator of all associated [`crate::annotations::OmimDisease`]s
    ///
    /// # Examples
    ///
    /// ```
    /// use hpo::{HpoTerm, Ontology};
    ///
    /// let ontology = Ontology::from_binary("tests/example.hpo").unwrap();
    ///
    /// let term = HpoTerm::try_new(&ontology, 11017u32.into()).unwrap();
    /// for disease in term.omim_diseases() {
    ///     println!("{}", disease.name());
    /// }
    /// ```
    pub fn omim_diseases(&self) -> OmimDiseaseIterator<'a> {
        OmimDiseaseIterator::new(self.omim_diseases, self.ontology)
    }

    /// Returns the set of all associated [`crate::annotations::OmimDisease`]s
    ///
    /// # Examples
    ///
    /// ```
    /// use hpo::{HpoTerm, Ontology};
    ///
    /// let ontology = Ontology::from_binary("tests/example.hpo").unwrap();
    ///
    /// let term = HpoTerm::try_new(&ontology, 1939u32.into()).unwrap();
    /// assert_eq!(term.omim_disease_ids().len(), 143);
    /// ```
    pub fn omim_disease_ids(&self) -> &OmimDiseases {
        self.omim_diseases
    }

    /// Returns the [`InformationContent`] of the term
    ///
    /// # Examples
    ///
    /// ```
    /// use hpo::{HpoTerm, Ontology};
    ///
    /// let ontology = Ontology::from_binary("tests/example.hpo").unwrap();
    ///
    /// let term = HpoTerm::try_new(&ontology, 1939u32.into()).unwrap();
    /// let ic = term.information_content();
    /// assert_eq!(ic.gene(), 0.6816717);
    /// assert_eq!(ic.omim_disease(), 3.4335358);
    /// ```
    pub fn information_content(&self) -> &InformationContent {
        self.information_content
    }

    /// Calculates the similarity of `self` and `other` using the provided `Similarity` algorithm
    ///
    /// # Examples
    ///
    /// ```
    /// use hpo::{HpoTerm, Ontology};
    /// use hpo::similarity::Builtins;
    /// use hpo::term::InformationContentKind;
    ///
    /// let ontology = Ontology::from_binary("tests/example.hpo").unwrap();
    ///
    /// let term1 = HpoTerm::try_new(&ontology, 11017u32.into()).unwrap();
    /// let term2 = HpoTerm::try_new(&ontology, 12639u32.into()).unwrap();
    ///
    /// let sim = term1.similarity_score(&term2, &Builtins::GraphIc(InformationContentKind::Omim));
    /// assert_eq!(sim, 0.2765914);
    /// ```
    pub fn similarity_score(&self, other: &HpoTerm, similarity: &impl Similarity) -> f32 {
        similarity.calculate(self, other)
    }

    /// Returns the distance (steps) from `self` to `other`, if `other` is a parent of `self`
    ///
    /// # Examples
    ///
    /// ```
    /// use hpo::{HpoTerm, Ontology};
    ///
    /// let ontology = Ontology::from_binary("tests/example.hpo").unwrap();
    ///
    /// let term1 = HpoTerm::try_new(&ontology, 25454u32.into()).unwrap();
    /// let term2 = HpoTerm::try_new(&ontology, 118u32.into()).unwrap();
    /// let term3 = HpoTerm::try_new(&ontology, 12639u32.into()).unwrap();
    ///
    /// assert_eq!(term1.distance_to_ancestor(&term2), Some(2));
    /// assert_eq!(term2.distance_to_ancestor(&term1), None);
    /// assert_eq!(term1.distance_to_ancestor(&term3), None);
    /// ```
    pub fn distance_to_ancestor(&self, other: &HpoTerm) -> Option<usize> {
        if self == other {
            return Some(0);
        }
        if self.parent_ids().contains(&other.id()) {
            return Some(1);
        }
        if !self.all_parent_ids().contains(&other.id()) {
            return None;
        }
        self.parents()
            .filter_map(|p| p.distance_to_ancestor(other))
            .min()
            .map(|c| c + 1)
    }

    /// Returns `true` if `self` is a child (direct or indirect) of `other`
    ///
    /// # Examples
    ///
    /// ```
    /// use hpo::{HpoTerm, Ontology};
    ///
    /// let ontology = Ontology::from_binary("tests/example.hpo").unwrap();
    ///
    /// let term1 = HpoTerm::try_new(&ontology, 25454u32.into()).unwrap();
    /// let term2 = HpoTerm::try_new(&ontology, 118u32.into()).unwrap();
    ///
    /// assert!(term1.child_of(&term2));
    /// assert!(!term2.child_of(&term1));
    /// ```
    pub fn child_of(&self, other: &HpoTerm) -> bool {
        self.all_parent_ids().contains(&other.id())
    }

    /// Returns `true` if `self` is a parent (direct or indirect) of `other`
    ///
    /// # Examples
    ///
    /// ```
    /// use hpo::{HpoTerm, Ontology};
    ///
    /// let ontology = Ontology::from_binary("tests/example.hpo").unwrap();
    ///
    /// let term1 = HpoTerm::try_new(&ontology, 25454u32.into()).unwrap();
    /// let term2 = HpoTerm::try_new(&ontology, 118u32.into()).unwrap();
    ///
    /// assert!(!term1.parent_of(&term2));
    /// assert!(term2.parent_of(&term1));
    /// ```
    pub fn parent_of(&self, other: &HpoTerm) -> bool {
        other.child_of(self)
    }

    /// Returns the shortest path to traverse from `self` ot `other`, if `other` is a parent of `self`
    ///
    /// # Examples
    ///
    /// ```
    /// use hpo::{HpoTerm, HpoTermId, Ontology};
    ///
    /// let ontology = Ontology::from_binary("tests/example.hpo").unwrap();
    ///
    /// let term1 = HpoTerm::try_new(&ontology, 25454u32.into()).unwrap();
    /// let term2 = HpoTerm::try_new(&ontology, 118u32.into()).unwrap();
    ///
    /// assert_eq!(
    ///     term1.path_to_ancestor(&term2).unwrap(),
    ///     vec![
    ///         HpoTermId::try_from("HP:0001939").unwrap(),
    ///         HpoTermId::try_from("HP:0000118").unwrap()
    ///     ]
    /// );
    /// assert_eq!(term2.path_to_ancestor(&term1), None);
    /// ```
    pub fn path_to_ancestor(&self, other: &HpoTerm) -> Option<Vec<HpoTermId>> {
        if self.id() == other.id() {
            return Some(vec![]);
        }
        if self.parent_ids().contains(&other.id()) {
            return Some(vec![other.id()]);
        }
        if !self.all_parent_ids().contains(&other.id()) {
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
    ///
    /// # Examples
    ///
    /// ```
    /// use hpo::{HpoTerm, Ontology};
    ///
    /// let ontology = Ontology::from_binary("tests/example.hpo").unwrap();
    ///
    /// let term1 = HpoTerm::try_new(&ontology, 25454u32.into()).unwrap();
    /// let term2 = HpoTerm::try_new(&ontology, 118u32.into()).unwrap();
    /// let term3 = HpoTerm::try_new(&ontology, 12639u32.into()).unwrap();
    ///
    /// assert_eq!(term1.distance_to_term(&term2), Some(2));
    /// assert_eq!(term2.distance_to_term(&term1), Some(2));
    /// assert_eq!(term1.distance_to_term(&term3), Some(4));
    /// ```
    pub fn distance_to_term(&self, other: &HpoTerm) -> Option<usize> {
        self.all_common_ancestors(other)
            .iter()
            .filter_map(|parent| {
                Some(self.distance_to_ancestor(&parent)? + other.distance_to_ancestor(&parent)?)
            })
            .min()
    }

    /// Returns `true` if the term is flagged as obsolete
    pub fn obsolete(&self) -> bool {
        self.obsolete
    }

    /// Returns the replacement term, if it exists
    ///
    /// Returns `None` otherwise
    ///
    /// If a term has multiple `replaced_by` annotations only one annotation
    /// is selected as replacement. Usually the last one, if parsing from an `obo` file.
    pub fn replaced_by(&self) -> Option<HpoTerm<'a>> {
        self.replaced_by
            .and_then(|term_id| self.ontology.hpo(term_id))
    }

    /// Returns the [`HpoTermId`] of an replacement term, if it exists
    ///
    /// Returns `None` otherwise
    ///
    /// If a term has multiple `replaced_by` annotations only one annotation
    /// is selected as replacement. Usually the last one, if parsing from an `obo` file.
    pub fn replacement_id(&self) -> Option<HpoTermId> {
        self.replaced_by
    }
}

impl PartialEq for HpoTerm<'_> {
    fn eq(&self, other: &Self) -> bool {
        self.id() == other.id()
    }
}

impl Eq for HpoTerm<'_> {}
