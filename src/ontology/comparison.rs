use std::collections::HashSet;
use std::fmt::Display;

use crate::annotations::{Gene, OmimDisease};
use crate::term::HpoGroup;
use crate::{HpoTerm, HpoTermId, Ontology};

#[derive(Debug)]
/// Compares the content of two Ontologies
///
/// This can be used when a new HPO masterdata release is available
/// to check what is changed between the previous and new one.
pub struct OntologyComparison<'a> {
    lhs: &'a Ontology,
    rhs: &'a Ontology,
}

impl<'a> Display for OntologyComparison<'a> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Terms\t{}\t{}\nGenes\t{}\t{}\nOmim Diseases\t{}\t{}",
            self.lhs.len(),
            self.rhs.len(),
            self.lhs.genes().count(),
            self.rhs.genes().count(),
            self.lhs.omim_diseases().count(),
            self.rhs.omim_diseases().count()
        )
    }
}

impl<'a> OntologyComparison<'a> {
    /// Constructs a new [`OntologyComparison`] from two [`Ontology`]
    /// The first argument, `lhs`, is considered the `old` or `base` Ontology,
    /// while the second argument, `rhs` is considered the `new` or `changed` one.
    pub fn new(lhs: &'a Ontology, rhs: &'a Ontology) -> Self {
        Self { lhs, rhs }
    }

    /// Returns all [`HpoTerm`]s that are exclusively in the `new` Ontology
    pub fn added_hpo_terms(&self) -> Vec<HpoTerm<'a>> {
        self.rhs
            .hpos()
            .filter(|term| self.lhs.hpo(term.id()).is_none())
            .collect()
    }

    /// Returns all [`HpoTerm`]s that are exclusively in the `old` Ontology
    pub fn removed_hpo_terms(&self) -> Vec<HpoTerm<'a>> {
        self.lhs
            .hpos()
            .filter(|term| self.rhs.hpo(term.id()).is_none())
            .collect()
    }

    /// Returns an [`HpoTermDelta`] struct for every HPO term that is different
    /// between the `old` and `new` Ontology.
    ///
    /// Differences are defined as either:
    /// - Changed name
    /// - Changed direct parents
    pub fn changed_hpo_terms(&self) -> Vec<HpoTermDelta> {
        self.lhs
            .hpos()
            .filter_map(|term| {
                if let Some(rhs) = self.rhs.hpo(term.id()) {
                    HpoTermDelta::new(term, rhs)
                } else {
                    None
                }
            })
            .collect()
    }

    /// Returns all [`Gene`]s that are exclusively in the `new` Ontology
    pub fn added_genes(&self) -> Vec<&Gene> {
        self.rhs
            .genes()
            .filter(|gene| self.lhs.gene(gene.id()).is_none())
            .collect()
    }

    /// Returns all [`Gene`]s that are exclusively in the `old` Ontology
    pub fn removed_genes(&self) -> Vec<&Gene> {
        self.lhs
            .genes()
            .filter(|gene| self.rhs.gene(gene.id()).is_none())
            .collect()
    }

    /// Returns an [`AnnotationDelta`] struct for every [`Gene`] that is different
    /// between the `old` and `new` Ontology.
    ///
    /// Differences are defined as either:
    /// - Changed name
    /// - Changed direct associated `HpoTerm`s
    pub fn changed_genes(&self) -> Vec<AnnotationDelta> {
        self.lhs
            .genes()
            .filter_map(|gene| {
                if let Some(rhs) = self.rhs.gene(gene.id()) {
                    AnnotationDelta::gene(gene, rhs)
                } else {
                    None
                }
            })
            .collect()
    }

    /// Returns all [`OmimDisease`]s that are exclusively in the `new` Ontology
    pub fn added_omim_diseases(&self) -> Vec<&OmimDisease> {
        self.rhs
            .omim_diseases()
            .filter(|disease| self.lhs.omim_disease(disease.id()).is_none())
            .collect()
    }

    /// Returns all [`OmimDisease`]s that are exclusively in the `old` Ontology
    pub fn removed_omim_diseases(&self) -> Vec<&OmimDisease> {
        self.lhs
            .omim_diseases()
            .filter(|disease| self.rhs.omim_disease(disease.id()).is_none())
            .collect()
    }

    /// Returns an [`AnnotationDelta`] struct for every [`OmimDisease`] that is different
    /// between the `old` and `new` Ontology.
    ///
    /// Differences are defined as either:
    /// - Changed name
    /// - Changed direct associated `HpoTerm`s
    pub fn changed_omim_diseases(&self) -> Vec<AnnotationDelta> {
        self.lhs
            .omim_diseases()
            .filter_map(|disease| {
                if let Some(rhs) = self.rhs.omim_disease(disease.id()) {
                    AnnotationDelta::disease(disease, rhs)
                } else {
                    None
                }
            })
            .collect()
    }
}

/// Differences between two [`HpoTerm`]s
pub struct HpoTermDelta {
    term_id: HpoTermId,
    changed_name: (String, String),
    added_parents: Vec<HpoTermId>,
    removed_parents: Vec<HpoTermId>,
}

impl HpoTermDelta {
    /// Constructs a new [`HpoTermDelta`] by comparing two [`HpoTerm`]s
    ///
    /// Returns `None` if both are identical
    pub fn new(lhs: HpoTerm, rhs: HpoTerm) -> Option<Self> {
        let changed_name = (lhs.name().to_string(), rhs.name().to_string());

        let lhs_parents: HashSet<HpoTermId> = lhs.parents().map(|t| t.id()).collect();
        let rhs_parents: HashSet<HpoTermId> = rhs.parents().map(|t| t.id()).collect();

        let removed_parents: Vec<HpoTermId> =
            lhs_parents.difference(&rhs_parents).copied().collect();
        let added_parents: Vec<HpoTermId> = rhs_parents.difference(&lhs_parents).copied().collect();

        if changed_name.0 != changed_name.1
            || !removed_parents.is_empty()
            || !added_parents.is_empty()
        {
            let term_id = lhs.id();
            Some(Self {
                term_id,
                changed_name,
                added_parents,
                removed_parents,
            })
        } else {
            None
        }
    }

    /// Returns all direct parent [`HpoTermId`]s of the `new` term that
    /// are not parents of the `old` term
    ///
    /// Returns `None` if no such terms exist
    pub fn added_parents(&self) -> Option<&Vec<HpoTermId>> {
        if self.added_parents.is_empty() {
            None
        } else {
            Some(&self.added_parents)
        }
    }

    /// Returns all direct parent [`HpoTermId`]s of the `old` term that
    /// are not parents of the `new` term
    ///
    /// Returns `None` if no such terms exist
    pub fn removed_parents(&self) -> Option<&Vec<HpoTermId>> {
        if self.removed_parents.is_empty() {
            None
        } else {
            Some(&self.removed_parents)
        }
    }

    /// Returns the `old` and `new` name if they are different
    ///
    /// Returns `None` if the name is unchanged
    pub fn changed_name(&self) -> Option<&(String, String)> {
        if self.changed_name.0 != self.changed_name.1 {
            Some(&self.changed_name)
        } else {
            None
        }
    }

    /// Returns the [`HpoTermId`] of the term
    pub fn id(&self) -> &HpoTermId {
        &self.term_id
    }
}

/// Differences between two [`Gene`]s or [`OmimDisease`]s
pub struct AnnotationDelta {
    id: String,
    names: (String, String),
    added_terms: Vec<HpoTermId>,
    removed_terms: Vec<HpoTermId>,
}

impl<'a> AnnotationDelta {
    /// Constructs a new [`AnnotationDelta`] by comparing two [`Gene`]s
    ///
    /// Returns `None` if both are identical
    pub fn gene(lhs: &Gene, rhs: &Gene) -> Option<Self> {
        let lhs_terms = lhs.hpo_terms();
        let rhs_terms = rhs.hpo_terms();

        let names = (lhs.name().to_string(), rhs.name().to_string());

        Self::delta(lhs_terms, rhs_terms, names, lhs.id().to_string())
    }

    /// Constructs a new [`AnnotationDelta`] by comparing two [`OmimDisease`]s
    ///
    /// Returns `None` if both are identical
    pub fn disease(lhs: &OmimDisease, rhs: &OmimDisease) -> Option<Self> {
        let lhs_terms = lhs.hpo_terms();
        let rhs_terms = rhs.hpo_terms();

        let names = (lhs.name().to_string(), rhs.name().to_string());

        Self::delta(lhs_terms, rhs_terms, names, lhs.id().to_string())
    }

    /// Calculates the difference between both items and constructs a new
    /// `AnnotationDelta` struct or returns None
    fn delta(
        lhs_terms: &HpoGroup,
        rhs_terms: &HpoGroup,
        names: (String, String),
        id: String,
    ) -> Option<Self> {
        let added_terms: Vec<HpoTermId> = rhs_terms
            .iter()
            .filter(|termid| !lhs_terms.contains(termid))
            .collect();

        let removed_terms: Vec<HpoTermId> = lhs_terms
            .iter()
            .filter(|termid| !rhs_terms.contains(termid))
            .collect();

        if !added_terms.is_empty() || !removed_terms.is_empty() || names.0 != names.1 {
            Some(Self {
                id,
                names,
                added_terms,
                removed_terms,
            })
        } else {
            None
        }
    }

    /// Returns all directly linked [`HpoTermId`]s of the `new` annotation that
    /// are not linked to the `old` anotation
    ///
    /// Returns `None` if no such terms exist
    pub fn added_terms(&'a self) -> Option<&Vec<HpoTermId>> {
        if self.added_terms.is_empty() {
            None
        } else {
            Some(&self.added_terms)
        }
    }

    /// Returns all directly linked [`HpoTermId`]s of the `old` annotation that
    /// are not linked to the `new` anotation
    ///
    /// Returns `None` if no such terms exist
    pub fn removed_terms(&'a self) -> Option<&Vec<HpoTermId>> {
        if self.removed_terms.is_empty() {
            None
        } else {
            Some(&self.removed_terms)
        }
    }

    /// Returns the `old` and `new` name if they are different
    ///
    /// Returns `None` if the name is unchanged
    pub fn changed_name(&self) -> Option<&(String, String)> {
        if self.names.0 != self.names.1 {
            Some(&self.names)
        } else {
            None
        }
    }

    /// Returns the `String`-formatted ID of the annotation
    pub fn id(&self) -> &str {
        &self.id
    }
}
