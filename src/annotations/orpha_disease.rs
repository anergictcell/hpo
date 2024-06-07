use std::collections::HashSet;
use std::fmt::{Debug, Display};
use std::hash::Hash;

use crate::annotations::disease::DiseaseIterator;
use crate::annotations::{AnnotationId, Disease};
use crate::term::HpoGroup;
use crate::HpoError;
use crate::HpoTermId;

/// A set of Orpha diseases
///
/// The set does not contain [`OrphaDisease`]s itself, but only
/// their [`OrphaDiseaseId`]s.
/// Currently implemented using [`HashSet`] but any other implementation
/// should work as well given that each [`OrphaDiseaseId`] must appear only once
/// and it provides an iterator of [`OrphaDiseaseId`]
pub type OrphaDiseases = HashSet<OrphaDiseaseId>;

/// A unique identifier for an [`OrphaDisease`]
///
/// This value can - in theory - represent any numerical unique value.
/// When using the default JAX provided masterdata, it represents
/// the actual Orpha MIM ID.
#[derive(Clone, Copy, Default, Debug, Hash, PartialEq, PartialOrd, Eq, Ord)]
pub struct OrphaDiseaseId {
    inner: u32,
}

impl AnnotationId for OrphaDiseaseId {
    /// Convert `self` to `u32`
    fn as_u32(&self) -> u32 {
        self.inner
    }
}

impl TryFrom<&str> for OrphaDiseaseId {
    type Error = HpoError;
    fn try_from(value: &str) -> Result<Self, Self::Error> {
        Ok(OrphaDiseaseId {
            inner: value.parse::<u32>()?,
        })
    }
}

impl From<u32> for OrphaDiseaseId {
    fn from(inner: u32) -> Self {
        OrphaDiseaseId { inner }
    }
}

impl Display for OrphaDiseaseId {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "ORPHA:{}", self.inner)
    }
}

/// A single ORPHA disease
///
/// A disease has a unique [`OrphaDiseaseId`] and a name and is
/// connected to a set of HPO terms
#[derive(Default, Debug, Clone)]
pub struct OrphaDisease {
    id: OrphaDiseaseId,
    name: String,
    hpos: HpoGroup,
}

impl Disease for OrphaDisease {
    type AnnoID = OrphaDiseaseId;

    /// Initializes a new Orpha disease
    ///
    /// This method should rarely, if ever, be used directly. The
    /// preferred way to create new genes is through [`Ontology::add_Orpha_disease`]
    /// to ensure that each disease exists only once.
    fn new(id: Self::AnnoID, name: &str) -> OrphaDisease {
        Self {
            name: name.to_string(),
            id,
            hpos: HpoGroup::default(),
        }
    }

    /// The unique [`OrphaDiseaseId`] of the disease, the Orpha MIM number
    fn id(&self) -> &Self::AnnoID {
        &self.id
    }

    /// The Orpha disease name
    fn name(&self) -> &str {
        &self.name
    }

    /// Connect another [HPO term](`crate::HpoTerm`) to the disease
    fn add_term<I: Into<HpoTermId>>(&mut self, term_id: I) -> bool {
        self.hpos.insert(term_id)
    }

    /// The set of connected HPO terms
    fn hpo_terms(&self) -> &HpoGroup {
        &self.hpos
    }
}

impl PartialEq for OrphaDisease {
    fn eq(&self, other: &OrphaDisease) -> bool {
        self.id == other.id
    }
}

impl Eq for OrphaDisease {}

impl TryFrom<&[u8]> for OrphaDisease {
    type Error = HpoError;
    /// Returns an [`OrphaDisease`] from a bytes vector
    ///
    /// The byte layout for this method is defined in
    /// [`Disease::as_bytes`]
    ///
    /// # Examples
    ///
    /// ```
    /// use hpo::annotations::{Disease, OrphaDisease, OrphaDiseaseId};
    ///
    /// let bytes = vec![
    ///     0u8, 0u8, 0u8, 22u8, // Total size of Blop
    ///     0u8, 0u8, 0u8, 123u8, // ID of the disease => 123
    ///     0u8, 0u8, 0u8, 6u8, // Length of name => 6
    ///     b'F', b'o', b'o', b'b', b'a', b'r', // Foobar
    ///     0u8, 0u8, 0u8, 0u8  // Number of associated HPO Terms => 0
    /// ];
    /// let disease = OrphaDisease::try_from(&bytes[..]).unwrap();
    ///
    /// assert_eq!(disease.name(), "Foobar");
    /// assert_eq!(disease.id(), &OrphaDiseaseId::from(123u32));
    /// ```
    fn try_from(bytes: &[u8]) -> Result<Self, Self::Error> {
        Self::from_bytes(bytes)
    }
}

impl Hash for OrphaDisease {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.id.hash(state);
    }
}

/// Iterates [`OrphaDisease`]
pub type OrphaDiseaseIterator<'a> = DiseaseIterator<'a, OrphaDiseaseId>;

impl<'a> std::iter::Iterator for DiseaseIterator<'a, OrphaDiseaseId> {
    type Item = &'a OrphaDisease;
    fn next(&mut self) -> Option<Self::Item> {
        self.diseases.next().map(|orpha_id| {
            self.ontology
                .orpha_disease(orpha_id)
                .expect("disease must exist in Ontology")
        })
    }
}
