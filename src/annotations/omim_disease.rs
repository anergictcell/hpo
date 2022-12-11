use std::collections::HashSet;
use std::fmt::Debug;
use std::fmt::Display;

use crate::term::HpoGroup;
use crate::HpoError;
use crate::HpoTermId;
use crate::Ontology;

/// A set of OMIM diseases
///
/// The set does not contain [`OmimDisease`]s itself, but only 
/// their [`OmimDiseaseId`]s.
/// Currently implemented using [`HashSet`] but any other implementation
/// should work as well given that each OmimDiseaseId must appear only once
/// and it provides an iterator of [`OmimDiseaseId`]
pub type OmimDiseases = HashSet<OmimDiseaseId>;


/// A unique identifier for an [`OmimDisease`]
///
/// This value can - in theory - represent any numerical unique value.
/// When using the default JAX provided masterdata, it represents
/// the actual OMIM MIM ID.
#[derive(Clone, Copy, Default, Debug, Hash, PartialEq, Eq)]
pub struct OmimDiseaseId {
    inner: usize,
}

impl TryFrom<&str> for OmimDiseaseId {
    type Error = HpoError;
    fn try_from(value: &str) -> Result<Self, Self::Error> {
        Ok(OmimDiseaseId {
            inner: value.parse::<usize>()?,
        })
    }
}

impl Display for OmimDiseaseId {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "OMIM:{}", self.inner)
    }
}

/// A single OMIM disease
///
/// A disease has a unique [`OmimDiseaseId`] and a name and is
/// connected to a set of HPO terms
pub struct OmimDisease {
    id: OmimDiseaseId,
    name: String,
    hpos: HpoGroup,
}

impl OmimDisease {
    /// Initializes a new OMIM disease
    ///
    /// This method should rarely, if ever, be used directly. The
    /// preferred way to create new genes is through [`Ontology::add_omim_disease`]
    /// to ensure that each disease exists only once.
    pub fn new(id: OmimDiseaseId, name: &str) -> OmimDisease {
        Self {
            name: name.to_string(),
            id,
            hpos: HpoGroup::default(),
        }
    }

    /// The unique [`OmimDiseaseId`] of the disease, the OMIM MIM number
    pub fn id(&self) -> &OmimDiseaseId {
        &self.id
    }

    /// The OMIM disease name
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Connect another [HPO term](`crate::HpoTerm`) to the disease 
    pub fn add_term(&mut self, term_id: HpoTermId) -> bool {
        self.hpos.insert(term_id)
    }

    /// The set of connected HPO terms
    pub fn hpo_terms(&self) -> &HpoGroup {
        &self.hpos
    }
}

impl PartialEq for OmimDisease {
    fn eq(&self, other: &OmimDisease) -> bool {
        self.id == other.id
    }
}
impl Eq for OmimDisease {}

/// [`OmimDisease`] Iterator
pub struct OmimDiseaseIterator<'a> {
    ontology: &'a Ontology,
    diseases: std::collections::hash_set::Iter<'a, OmimDiseaseId>,
}

impl<'a> OmimDiseaseIterator<'a> {
    /// Initialize a new [`OmimDiseaseIterator`]
    ///
    /// This method requires the [`Ontology`] as a parameter since
    /// the actual [`OmimDisease`] entities are stored in it and
    /// not in [`OmimDiseases`] itself
    pub fn new(diseases: &'a OmimDiseases, ontology: &'a Ontology) -> Self {
        OmimDiseaseIterator {
            diseases: diseases.iter(),
            ontology,
        }
    }
}

impl<'a> std::iter::Iterator for OmimDiseaseIterator<'a> {
    type Item = &'a OmimDisease;
    fn next(&mut self) -> Option<Self::Item> {
        match self.diseases.next() {
            Some(omim_id) => Some(self.ontology.omim_disease(omim_id).unwrap()),
            None => None,
        }
    }
}

impl Debug for OmimDiseaseIterator<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "OmimDiseaseIterator")
    }
}
