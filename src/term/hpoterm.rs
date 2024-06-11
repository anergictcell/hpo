use crate::annotations::GeneIterator;
use crate::annotations::Genes;
use crate::annotations::OmimDiseaseIterator;
use crate::annotations::OmimDiseases;
use crate::annotations::OrphaDiseaseIterator;
use crate::annotations::OrphaDiseases;
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
    orpha_diseases: &'a OrphaDiseases,
    information_content: &'a InformationContent,
    obsolete: bool,
    replaced_by: Option<HpoTermId>,
    ontology: &'a Ontology,
}

impl<'a> HpoTerm<'a> {
    /// Constructs a new [`HpoTerm`]
    ///
    /// Clients should rarely, if ever, use this function. Instead, use the
    /// [`Ontology::hpo`] method to get an [`HpoTerm`]
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
    /// let term = HpoTerm::try_new(&ontology, 118u32);
    /// assert!(term.is_ok());
    ///
    /// let non_existing_term = HpoTerm::try_new(&ontology, 666666666u32);
    /// assert!(non_existing_term.is_err());
    ///
    /// // better retrieve terms from `Ontology`:
    ///
    /// let term = ontology.hpo(118u32);
    /// assert!(term.is_some());
    /// ```
    pub fn try_new<I: Into<HpoTermId>>(ontology: &'a Ontology, term: I) -> HpoResult<HpoTerm<'a>> {
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
            orpha_diseases: term.orpha_diseases(),
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
    /// let term = ontology.hpo(118u32).unwrap();
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
    /// let term = ontology.hpo(118u32).unwrap();
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
    /// let term = ontology.hpo(11017u32).unwrap();
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
    /// let term = ontology.hpo(11017u32).unwrap();
    /// assert_eq!(term.parents().count(), 1);
    /// for parent in term.parents() {
    ///    println!("{}", parent.name());
    /// }
    /// ```
    pub fn parents(&self) -> Iter<'a> {
        Iter::new(self.parents.iter(), self.ontology)
    }

    /// Returns the [`HpoTermId`]s of all direct and indirect parents
    ///
    /// # Examples
    ///
    /// ```
    /// use hpo::{HpoTerm, Ontology};
    ///
    /// let ontology = Ontology::from_binary("tests/example.hpo").unwrap();
    ///
    /// let term = ontology.hpo(11017u32).unwrap();
    /// assert_eq!(term.all_parent_ids().len(), 3);
    /// ```
    pub fn all_parent_ids(&self) -> &HpoGroup {
        self.all_parents
    }

    /// Returns an iterator of the direct and indirect parents of the term
    ///
    /// # Examples
    ///
    /// ```
    /// use hpo::{HpoTerm, Ontology};
    ///
    /// let ontology = Ontology::from_binary("tests/example.hpo").unwrap();
    ///
    /// let term = ontology.hpo(11017u32).unwrap();
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
    /// let term = ontology.hpo(1939u32).unwrap();
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
    /// let term = ontology.hpo(1939u32).unwrap();
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
    /// let term1 = ontology.hpo(11017u32).unwrap();
    /// let term2 = ontology.hpo(25454u32).unwrap();
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
    /// let term1 = ontology.hpo(11017u32).unwrap();
    /// let term2 = ontology.hpo(1939u32).unwrap();
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
    /// let term1 = ontology.hpo(11017u32).unwrap();
    /// let term2 = ontology.hpo(25454u32).unwrap();
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
    /// let term1 = ontology.hpo(11017u32).unwrap();
    /// let term2 = ontology.hpo(1939u32).unwrap();
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
    /// let term1 = ontology.hpo(11017u32).unwrap();
    /// let term2 = ontology.hpo(12639u32).unwrap();
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
    /// let term1 = ontology.hpo(11017u32).unwrap();
    /// let term2 = ontology.hpo(1939u32).unwrap();
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
    /// let term1 = ontology.hpo(11017u32).unwrap();
    /// let term2 = ontology.hpo(12639u32).unwrap();
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
    /// let term1 = ontology.hpo(11017u32).unwrap();
    /// let term2 = ontology.hpo(1939u32).unwrap();
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
    /// let term1 = ontology.hpo(11017u32).unwrap();
    /// let term2 = ontology.hpo(25454u32).unwrap();
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
    /// let term1 = ontology.hpo(11017u32).unwrap();
    /// let term2 = ontology.hpo(1939u32).unwrap();
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
    /// let term1 = ontology.hpo(11017u32).unwrap();
    /// let term2 = ontology.hpo(25454u32).unwrap();
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
    /// let term1 = ontology.hpo(11017u32).unwrap();
    /// let term2 = ontology.hpo(1939u32).unwrap();
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
    /// let term1 = ontology.hpo(11017u32).unwrap();
    /// let term2 = ontology.hpo(12639u32).unwrap();
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
    /// let term1 = ontology.hpo(11017u32).unwrap();
    /// let term2 = ontology.hpo(1939u32).unwrap();
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
    /// let term1 = ontology.hpo(11017u32).unwrap();
    /// let term2 = ontology.hpo(12639u32).unwrap();
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
    /// let term1 = ontology.hpo(11017u32).unwrap();
    /// let term2 = ontology.hpo(1939u32).unwrap();
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
    /// let term = ontology.hpo(11017u32).unwrap();
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
    /// let term = ontology.hpo(12638u32).unwrap();
    /// assert_eq!(term.gene_ids().len(), 41);
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
    /// use hpo::annotations::Disease;
    ///
    /// let ontology = Ontology::from_binary("tests/example.hpo").unwrap();
    ///
    /// let term = ontology.hpo(11017u32).unwrap();
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
    /// let term = ontology.hpo(1939u32).unwrap();
    /// assert_eq!(term.omim_disease_ids().len(), 93);
    /// ```
    pub fn omim_disease_ids(&self) -> &OmimDiseases {
        self.omim_diseases
    }

    /// Returns an iterator of all associated [`crate::annotations::OrphaDisease`]s
    ///
    /// # Examples
    ///
    /// ```
    /// use hpo::{HpoTerm, Ontology};
    /// use hpo::annotations::Disease;
    ///
    /// let ontology = Ontology::from_binary("tests/example.hpo").unwrap();
    ///
    /// let term = ontology.hpo(11017u32).unwrap();
    /// for disease in term.orpha_diseases() {
    ///     println!("{}", disease.name());
    /// }
    /// ```
    pub fn orpha_diseases(&self) -> OrphaDiseaseIterator<'a> {
        OrphaDiseaseIterator::new(self.orpha_diseases, self.ontology)
    }

    /// Returns the set of all associated [`crate::annotations::OrphaDisease`]s
    ///
    /// # Examples
    ///
    /// ```
    /// use hpo::{HpoTerm, Ontology};
    /// use hpo::annotations::Disease;
    ///
    /// let ontology = Ontology::from_binary("tests/example.hpo").unwrap();
    ///
    /// let term = ontology.hpo(1939u32).unwrap();
    /// assert_eq!(term.orpha_disease_ids().len(), 32);
    /// ```
    pub fn orpha_disease_ids(&self) -> &OrphaDiseases {
        self.orpha_diseases
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
    /// let term = ontology.hpo(1939u32).unwrap();
    /// let ic = term.information_content();
    /// assert_eq!(ic.gene(), 1.9442855);
    /// assert_eq!(ic.omim_disease(), 0.4578331);
    /// assert_eq!(ic.orpha_disease(), 2.2994552);
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
    /// let term1 = ontology.hpo(12638u32).unwrap(); // Abnormal nervous system physiology
    /// let term2 = ontology.hpo(100547u32).unwrap(); // Abnormal forebrain morphology
    ///
    /// let sim = term1.similarity_score(&term2, &Builtins::GraphIc(InformationContentKind::Omim));
    /// assert_eq!(sim, 0.112043366);
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
    /// let term1 = ontology.hpo(25454u32).unwrap();
    /// let term2 = ontology.hpo(118u32).unwrap();
    /// let term3 = ontology.hpo(12639u32).unwrap();
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
    /// let term1 = ontology.hpo(25454u32).unwrap();
    /// let term2 = ontology.hpo(118u32).unwrap();
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
    /// let term1 = ontology.hpo(25454u32).unwrap();
    /// let term2 = ontology.hpo(118u32).unwrap();
    ///
    /// assert!(!term1.parent_of(&term2));
    /// assert!(term2.parent_of(&term1));
    /// ```
    pub fn parent_of(&self, other: &HpoTerm) -> bool {
        other.child_of(self)
    }

    /// Returns the shortest path to traverse from `self` to `other`, if `other` is a parent of `self`
    ///
    /// # Examples
    ///
    /// ```
    /// use hpo::{HpoTerm, HpoTermId, Ontology};
    ///
    /// let ontology = Ontology::from_binary("tests/example.hpo").unwrap();
    ///
    /// let term1 = ontology.hpo(25454u32).unwrap();
    /// let term2 = ontology.hpo(118u32).unwrap();
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
    /// let term1 = ontology.hpo(25454u32).unwrap();
    /// let term2 = ontology.hpo(118u32).unwrap();
    /// let term3 = ontology.hpo(12639u32).unwrap();
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

    /// Returns the shortest path to traverse from `self` to `other`
    ///
    /// This method is not optimized for performance
    ///
    /// # Examples
    ///
    /// ```
    /// use hpo::{HpoTerm, Ontology, HpoTermId};
    ///
    /// let ontology = Ontology::from_binary("tests/example.hpo").unwrap();
    ///
    /// let term1 = ontology.hpo(25454u32).unwrap();
    /// let term3 = ontology.hpo(12639u32).unwrap();
    ///
    /// assert_eq!(
    ///     term1.path_to_term(&term3).unwrap(),
    ///     vec![
    ///         HpoTermId::try_from("HP:0001939").unwrap(),
    ///         HpoTermId::try_from("HP:0000118").unwrap(),
    ///         HpoTermId::try_from("HP:0000707").unwrap(),
    ///         HpoTermId::try_from("HP:0012639").unwrap()
    ///     ]
    /// );
    /// ```
    /// # Panics
    /// TODO    
    pub fn path_to_term(&self, other: &HpoTerm) -> Option<Vec<HpoTermId>> {
        if other.parent_of(self) {
            return self.path_to_ancestor(other);
        }
        if self.parent_of(other) {
            return other.path_to_ancestor(self).map(|terms| {
                terms
                    .iter()
                    .rev()
                    .skip(1)
                    .chain(std::iter::once(&other.id()))
                    .copied()
                    .collect()
            });
        }

        self.all_common_ancestors(other)
            .iter()
            .map(|ancestor| {
                (
                    ancestor,
                    self.distance_to_ancestor(&ancestor)
                        .expect("self must have a path to its ancestor")
                        + other
                            .distance_to_ancestor(&ancestor)
                            .expect("other must have a path to its ancestor"),
                )
            })
            .min_by_key(|tuple| tuple.1)
            .map(|min| {
                self.path_to_ancestor(&min.0)
                    .expect("self must have a path to its ancestor")
                    .iter()
                    .chain(
                        other
                            .path_to_ancestor(&min.0)
                            .expect("other must have a path to its ancestor")
                            .iter()
                            .rev()
                            .skip(1),
                    )
                    .chain(std::iter::once(&other.id()))
                    .copied()
                    .collect()
            })
    }

    /// Returns `true` if the term is flagged as obsolete
    pub fn is_obsolete(&self) -> bool {
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

    /// Returns `true` if the term is a descendent of a modifier root term
    ///
    /// # Examples
    ///
    /// ```
    /// use hpo::{HpoTerm, Ontology};
    ///
    /// let ontology = Ontology::from_binary("tests/example.hpo").unwrap();
    ///
    /// let mendelian_inheritance = ontology.hpo(34345u32).unwrap();
    /// let adult_onset = ontology.hpo(3581u32).unwrap();
    /// let abnormal_forebrain_morphology = ontology.hpo(100547u32).unwrap();
    ///
    /// assert!(mendelian_inheritance.is_modifier());
    /// assert!(adult_onset.is_modifier());
    /// assert!(!abnormal_forebrain_morphology.is_modifier());
    /// ```
    pub fn is_modifier(&self) -> bool {
        self.ontology
            .modifier()
            .iter()
            .any(|modifier_root| (self.all_parent_ids() | self.id()).contains(&modifier_root))
    }

    /// Returns the category of the term, or `None` if uncategorized
    ///
    /// Categories are defined in [`Ontology::set_default_categories()`]
    ///
    /// The default implementation ensures that each term is categorized
    /// in exactly one category, but this can be modified on the Ontology level.
    /// If a term might be part of multiple categories, one is (semi-)randomly
    /// selected.
    ///
    /// # Examples
    ///
    /// ```
    /// use hpo::{HpoTerm, Ontology, HpoTermId};
    ///
    /// let ontology = Ontology::from_binary("tests/example.hpo").unwrap();
    ///
    /// let mendelian_inheritance = ontology.hpo(HpoTermId::from_u32(34345)).unwrap();
    /// let adult_onset = ontology.hpo(HpoTermId::from_u32(3581)).unwrap();
    /// let abnormal_hypothalamus_physiology = ontology.hpo(HpoTermId::from_u32(12285)).unwrap();
    ///
    /// assert_eq!(mendelian_inheritance.categories(), vec![HpoTermId::from_u32(5)]);
    /// assert_eq!(adult_onset.categories(), vec![HpoTermId::from_u32(12823)]);
    /// assert_eq!(abnormal_hypothalamus_physiology.categories(), vec![HpoTermId::from_u32(707), HpoTermId::from_u32(818)]);
    /// ```
    pub fn categories(&self) -> Vec<HpoTermId> {
        self.ontology
            .categories()
            .iter()
            .filter(|cat| (self.all_parent_ids() | self.id()).contains(cat))
            .collect()
    }
}

impl PartialEq for HpoTerm<'_> {
    fn eq(&self, other: &Self) -> bool {
        self.id() == other.id()
    }
}

impl Eq for HpoTerm<'_> {}

#[cfg(test)]
mod test_categories {
    use super::*;
    use crate::Ontology;

    #[test]
    fn test_categories() {
        let ontology = Ontology::from_binary("tests/example.hpo").unwrap();

        let mendelian_inheritance = ontology.hpo(HpoTermId::from_u32(34345)).unwrap();
        let adult_onset = ontology.hpo(HpoTermId::from_u32(3581)).unwrap();
        let abnormal_hypothalamus_physiology = ontology.hpo(HpoTermId::from_u32(12285)).unwrap();

        assert_eq!(
            mendelian_inheritance.categories(),
            vec![HpoTermId::from_u32(5)]
        );
        assert_eq!(adult_onset.categories(), vec![HpoTermId::from_u32(12823)]);
        assert_eq!(
            abnormal_hypothalamus_physiology.categories(),
            vec![HpoTermId::from_u32(707), HpoTermId::from_u32(818)]
        );
    }

    #[test]
    fn test_categories_is_self() {
        let ontology = Ontology::from_binary("tests/example.hpo").unwrap();

        let inheritance = ontology.hpo(HpoTermId::from_u32(5)).unwrap();
        assert_eq!(inheritance.categories(), vec![HpoTermId::from_u32(5)]);

        let nervous_system = ontology.hpo(HpoTermId::from_u32(707)).unwrap();
        assert_eq!(nervous_system.categories(), vec![HpoTermId::from_u32(707)]);
    }

    #[test]
    fn test_modifier() {
        let ontology = Ontology::from_binary("tests/example.hpo").unwrap();

        let mendelian_inheritance = ontology.hpo(34345u32).unwrap();
        let adult_onset = ontology.hpo(3581u32).unwrap();
        let abnormal_forebrain_morphology = ontology.hpo(100_547u32).unwrap();

        assert!(mendelian_inheritance.is_modifier());
        assert!(adult_onset.is_modifier());
        assert!(!abnormal_forebrain_morphology.is_modifier());
    }

    #[test]
    fn test_modifier_is_self() {
        let ontology = Ontology::from_binary("tests/example.hpo").unwrap();

        let inheritance = ontology.hpo(5u32).unwrap();
        let nervous_system = ontology.hpo(707u32).unwrap();

        assert!(inheritance.is_modifier());
        assert!(!nervous_system.is_modifier());
    }
}

#[cfg(test)]
mod test_path_to_term {
    use super::*;

    #[test]
    fn normal() {
        let ontology = Ontology::from_binary("tests/example.hpo").unwrap();
        let term1 = ontology.hpo(25454u32).unwrap();
        let term2 = ontology.hpo(12639u32).unwrap();

        assert_eq!(
            term1.path_to_term(&term2).unwrap(),
            vec![
                HpoTermId::try_from("HP:0001939").unwrap(),
                HpoTermId::try_from("HP:0000118").unwrap(),
                HpoTermId::try_from("HP:0000707").unwrap(),
                HpoTermId::try_from("HP:0012639").unwrap()
            ]
        );

        assert_eq!(
            term2.path_to_term(&term1).unwrap(),
            vec![
                HpoTermId::try_from("HP:0000707").unwrap(),
                HpoTermId::try_from("HP:0000118").unwrap(),
                HpoTermId::try_from("HP:0001939").unwrap(),
                HpoTermId::try_from("HP:0025454").unwrap(),
            ]
        );
    }

    #[test]
    fn same_term() {
        let ontology = Ontology::from_binary("tests/example.hpo").unwrap();
        let term1 = ontology.hpo(12639u32).unwrap();
        let term2 = ontology.hpo(12639u32).unwrap();

        assert_eq!(
            term1.path_to_term(&term2).unwrap(),
            vec![HpoTermId::try_from("HP:0012639").unwrap()]
        );
    }

    #[test]
    fn parent_term() {
        let ontology = Ontology::from_binary("tests/example.hpo").unwrap();
        let term1 = ontology.hpo(12639u32).unwrap();
        let term2 = ontology.hpo(707u32).unwrap();

        assert_eq!(
            term1.path_to_term(&term2).unwrap(),
            vec![HpoTermId::try_from("HP:0000707").unwrap()]
        );
    }

    #[test]
    fn child_term() {
        let ontology = Ontology::from_binary("tests/example.hpo").unwrap();
        let term1 = ontology.hpo(707u32).unwrap();
        let term2 = ontology.hpo(12639u32).unwrap();

        assert_eq!(
            term1.path_to_term(&term2).unwrap(),
            vec![HpoTermId::try_from("HP:0012639").unwrap()]
        );
    }
}
