use std::collections::hash_map::Values;
use std::collections::HashSet;
use std::fmt::{Debug, Display};
use std::hash::Hash;

use crate::annotations::disease::DiseaseIterator;
use crate::annotations::{AnnotationId, Disease};
use crate::term::HpoGroup;
use crate::{HpoError, HpoSet, HpoTermId, Ontology};

/// A set of OMIM diseases
///
/// The set does not contain [`OmimDisease`]s itself, but only
/// their [`OmimDiseaseId`]s.
/// Currently implemented using [`HashSet`] but any other implementation
/// should work as well given that each [`OmimDiseaseId`] must appear only once
/// and it provides an iterator of [`OmimDiseaseId`]
pub type OmimDiseases = HashSet<OmimDiseaseId>;

/// A unique identifier for an [`OmimDisease`]
///
/// This value can - in theory - represent any numerical unique value.
/// When using the default JAX provided masterdata, it represents
/// the actual OMIM MIM ID.
#[derive(Clone, Copy, Default, Debug, Hash, PartialEq, PartialOrd, Eq, Ord)]
pub struct OmimDiseaseId {
    inner: u32,
}

impl AnnotationId for OmimDiseaseId {
    /// Convert `self` to `u32`
    fn as_u32(&self) -> u32 {
        self.inner
    }
}

impl TryFrom<&str> for OmimDiseaseId {
    type Error = HpoError;
    fn try_from(value: &str) -> Result<Self, Self::Error> {
        Ok(OmimDiseaseId {
            inner: value.parse::<u32>()?,
        })
    }
}

impl From<u32> for OmimDiseaseId {
    fn from(inner: u32) -> Self {
        OmimDiseaseId { inner }
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
#[derive(Default, Debug, Clone)]
pub struct OmimDisease {
    id: OmimDiseaseId,
    name: String,
    hpos: HpoGroup,
}

impl Disease for OmimDisease {
    type AnnoID = OmimDiseaseId;

    /// Initializes a new OMIM disease
    ///
    /// This method should rarely, if ever, be used directly. The
    /// preferred way to create new genes is through [`Builder::annotate_omim_disease`](`crate::builder::Builder::annotate_omim_disease`)
    /// to ensure that each disease exists only once.
    fn new(id: Self::AnnoID, name: &str) -> OmimDisease {
        Self {
            name: name.to_string(),
            id,
            hpos: HpoGroup::default(),
        }
    }

    /// The unique [`OmimDiseaseId`] of the disease, the OMIM MIM number
    fn id(&self) -> &Self::AnnoID {
        &self.id
    }

    /// The OMIM disease name
    fn name(&self) -> &str {
        &self.name
    }

    /// The set of connected HPO terms
    fn hpo_terms(&self) -> &HpoGroup {
        &self.hpos
    }

    /// Returns a binary representation of the `OmimDisease`
    ///
    /// The binary layout is defined as:
    ///
    /// | Byte offset | Number of bytes | Description |
    /// | --- | --- | --- |
    /// | 0 | 4 | The total length of the binary data blob as big-endian `u32` |
    /// | 4 | 4 | The `OmimDiseaseId` as big-endian `u32` |
    /// | 8 | 4 | The length of the `OmimDisease` Name as big-endian `u32` |
    /// | 12 | n | The `OmimDisease` name as u8 vector |
    /// | 12 + n | 4 | The number of associated HPO terms as big-endian `u32` |
    /// | 16 + n | x * 4 | The [`HpoTermId`]s of the associated terms, each encoded as big-endian `u32` |
    ///
    /// # Examples
    ///
    /// ```
    /// use hpo::annotations::{Disease, OmimDisease};
    ///
    /// let mut disease = OmimDisease::new(123.into(), "FooBar");
    /// let bytes = disease.as_bytes();
    ///
    /// assert_eq!(bytes.len(), 4 + 4 + 4 + 6 + 4);
    /// assert_eq!(bytes[4..8], [0u8, 0u8, 0u8, 123u8]); // ID of disease => 123
    /// assert_eq!(bytes[8..12], [0u8, 0u8, 0u8, 6u8]); // Length of Name => 6
    /// ```
    fn as_bytes(&self) -> Vec<u8> {
        fn usize_to_u32(n: usize) -> u32 {
            n.try_into().expect("unable to convert {n} to u32")
        }
        let name = self.name().as_bytes();
        let name_length = name.len();
        let size = 4 + 4 + 4 + name_length + 4 + self.hpos.len() * 4;

        let mut res = Vec::new();

        // 4 bytes for total length
        res.append(&mut usize_to_u32(size).to_be_bytes().to_vec());

        // 4 bytes for OMIM Disease-ID
        res.append(&mut self.id.to_be_bytes().to_vec());

        // 4 bytes for Length of OMIM Disease Name
        res.append(&mut usize_to_u32(name_length).to_be_bytes().to_vec());

        // OMIM Disease name (n bytes)
        for c in name {
            res.push(*c);
        }

        // 4 bytes for number of HPO terms
        res.append(&mut usize_to_u32(self.hpos.len()).to_be_bytes().to_vec());

        // HPO terms
        res.append(&mut self.hpos.as_bytes());

        res
    }

    /// Returns an [`HpoSet`] from the `OmimDisease`
    fn to_hpo_set<'a>(&self, ontology: &'a Ontology) -> HpoSet<'a> {
        HpoSet::new(ontology, self.hpos.clone())
    }

    /// Connect another [HPO term](`crate::HpoTerm`) to the disease
    ///
    /// # Note
    ///
    /// This method does **not** add the [`OmimDisease`] to the [HPO term](`crate::HpoTerm`).
    /// Clients should not use this method, unless they are creating their own Ontology.
    fn add_term<I: Into<HpoTermId>>(&mut self, term_id: I) -> bool {
        self.hpos.insert(term_id)
    }
}

impl PartialEq for OmimDisease {
    fn eq(&self, other: &OmimDisease) -> bool {
        self.id == other.id
    }
}

impl Eq for OmimDisease {}

impl TryFrom<&[u8]> for OmimDisease {
    type Error = HpoError;
    /// Returns an [`OmimDisease`] from a bytes vector
    ///
    /// The byte layout for this method is defined in
    /// [`Disease::as_bytes`]
    ///
    /// # Examples
    ///
    /// ```
    /// use hpo::annotations::{Disease, OmimDisease, OmimDiseaseId};
    ///
    /// let bytes = vec![
    ///     0u8, 0u8, 0u8, 22u8, // Total size of Blop
    ///     0u8, 0u8, 0u8, 123u8, // ID of the disease => 123
    ///     0u8, 0u8, 0u8, 6u8, // Length of name => 6
    ///     b'F', b'o', b'o', b'b', b'a', b'r', // Foobar
    ///     0u8, 0u8, 0u8, 0u8  // Number of associated HPO Terms => 0
    /// ];
    /// let disease = OmimDisease::try_from(&bytes[..]).unwrap();
    ///
    /// assert_eq!(disease.name(), "Foobar");
    /// assert_eq!(disease.id(), &OmimDiseaseId::from(123u32));
    /// ```
    fn try_from(bytes: &[u8]) -> Result<Self, Self::Error> {
        Self::from_bytes(bytes)
    }
}

impl Hash for OmimDisease {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.id.hash(state);
    }
}

/// Iterates [`OmimDisease`]
pub type OmimDiseaseIterator<'a> = DiseaseIterator<'a, OmimDiseaseId>;

impl<'a> std::iter::Iterator for DiseaseIterator<'a, OmimDiseaseId> {
    type Item = &'a OmimDisease;
    fn next(&mut self) -> Option<Self::Item> {
        self.diseases.next().map(|omim_id| {
            self.ontology
                .omim_disease(omim_id)
                .expect("disease must exist in Ontology")
        })
    }
}

/// Iterates [`OmimDisease`] that match the query string
///
/// This struct is returned by [`crate::Ontology::omim_diseases_by_name`]
pub struct OmimDiseaseFilter<'a> {
    iter: Values<'a, OmimDiseaseId, OmimDisease>,
    query: &'a str,
}

impl<'a> OmimDiseaseFilter<'a> {
    pub(crate) fn new(iter: Values<'a, OmimDiseaseId, OmimDisease>, query: &'a str) -> Self {
        OmimDiseaseFilter { iter, query }
    }
}

impl<'a> Iterator for OmimDiseaseFilter<'a> {
    type Item = &'a OmimDisease;

    fn next(&mut self) -> Option<Self::Item> {
        self.iter
            .by_ref()
            .find(|&item| item.name().contains(self.query))
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn disease_to_binary() {
        let mut disease = OmimDisease::new(123u32.into(), "FooBar");
        disease.add_term(66u32);
        disease.add_term(77u32);

        let bin = disease.as_bytes();

        assert_eq!(bin.len(), 4 + 4 + 4 + 6 + 4 + 8);
    }

    #[test]
    fn disease_to_and_from_binary() {
        let mut disease = OmimDisease::new(123u32.into(), "FooBar");
        disease.add_term(66u32);
        disease.add_term(77u32);

        let bin = disease.as_bytes();

        let disease2 = OmimDisease::try_from(&bin[..]).expect("Can't build Disease");

        assert_eq!(disease.name(), disease2.name());
        assert_eq!(disease.id(), disease2.id());
        assert_eq!(disease.hpo_terms().len(), disease2.hpo_terms().len());
        for (a, b) in disease.hpo_terms().iter().zip(disease2.hpo_terms().iter()) {
            assert_eq!(a, b);
        }
    }

    #[test]
    fn disease_with_non_utf8_name() {
        let mut disease = OmimDisease::new(123u32.into(), "Foo😀Bar");
        disease.add_term(66u32);
        disease.add_term(77u32);

        let bin = disease.as_bytes();

        // the smiley uses 4 bytes, so the name length is 10
        assert_eq!(bin.len(), 4 + 4 + 4 + 10 + 4 + 8);

        let disease2 = OmimDisease::try_from(&bin[..]).expect("Can't build Disease");

        assert_eq!(disease.name(), disease2.name());
        assert_eq!(disease.id(), disease2.id());
        assert_eq!(disease.hpo_terms().len(), disease2.hpo_terms().len());
        for (a, b) in disease.hpo_terms().iter().zip(disease2.hpo_terms().iter()) {
            assert_eq!(a, b);
        }
    }

    #[test]
    fn disease_with_long_name() {
        let name = "Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat. Duis aute irure dolor in reprehenderit in voluptate velit esse cillum dolore eu fugiat nulla pariatur. Excepteur sint occaecat cupidatat non proident, sunt in culpa qui officia deserunt mollit anim id est laborum.";
        let mut disease = OmimDisease::new(123u32.into(), name);
        disease.add_term(66u32);
        disease.add_term(77u32);

        let bin = disease.as_bytes();

        let disease2 = OmimDisease::try_from(&bin[..]).expect("Can't build Disease");

        assert_eq!(disease.name(), disease2.name());
        assert_eq!(disease.id(), disease2.id());
        assert_eq!(disease.hpo_terms().len(), disease2.hpo_terms().len());
        for (a, b) in disease.hpo_terms().iter().zip(disease2.hpo_terms().iter()) {
            assert_eq!(a, b);
        }
    }

    #[test]
    fn disease_with_wrong_length() {
        let mut disease = OmimDisease::new(123u32.into(), "foobar");
        disease.add_term(66u32);
        disease.add_term(77u32);

        let mut bin = disease.as_bytes();

        assert!(OmimDisease::try_from(&bin[..15]).is_err());
        assert!(OmimDisease::try_from(&bin[..21]).is_err());
        assert!(OmimDisease::try_from(&bin[..29]).is_err());
        assert!(OmimDisease::try_from(&bin[..30]).is_ok());

        bin.push(1);
        assert!(OmimDisease::try_from(&bin[..30]).is_ok());
        assert!(OmimDisease::try_from(&bin[..31]).is_err());
    }
}
