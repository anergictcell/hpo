use std::collections::HashSet;
use std::fmt::Debug;
use std::hash::Hash;
use tracing::error;

use crate::annotations::AnnotationId;
use crate::term::HpoGroup;
use crate::u32_from_bytes;
use crate::{HpoError, HpoSet, HpoTermId, Ontology};

/// Defines common methods that all types of diseases have
pub trait Disease: PartialEq + Eq + Hash + Sized {
    /// The type of the corresponding ID
    /// e.g. `OmimDisease` has `OmimDiseaseId`
    type AnnoID: AnnotationId;

    /// Initializes a new disease
    fn new(id: Self::AnnoID, name: &str) -> Self;

    /// Returns the ID of the disease
    fn id(&self) -> &Self::AnnoID;

    /// Returns the name of the disease
    fn name(&self) -> &str;

    /// Connect another [HPO term](`crate::HpoTerm`) to the disease
    fn add_term<I: Into<HpoTermId>>(&mut self, term_id: I) -> bool;

    /// Returns a reference to the associated [group of HPO terms](`crate::term::HpoGroup`)
    fn hpo_terms(&self) -> &HpoGroup;

    /// Creates a new `crate::set::HpoSet`
    fn to_hpo_set<'a>(&self, ontology: &'a Ontology) -> HpoSet<'a> {
        HpoSet::new(ontology, self.hpo_terms().clone())
    }

    /// Returns a binary representation of the `Disease`
    ///
    /// The binary layout is defined as:
    ///
    /// | Byte offset | Number of bytes | Description |
    /// | --- | --- | --- |
    /// | 0 | 4 | The total length of the binary data blob as big-endian `u32` |
    /// | 4 | 4 | The `OmimDiseaseId` as big-endian `u32` |
    /// | 8 | 4 | The length of the `Disease` Name as big-endian `u32` |
    /// | 12 | n | The Disease name as u8 vector |
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
        let size = 4 + 4 + 4 + name_length + 4 + self.hpo_terms().len() * 4;

        let mut res = Vec::new();

        // 4 bytes for total length
        res.append(&mut usize_to_u32(size).to_be_bytes().to_vec());

        // 4 bytes for Disease-ID
        res.append(&mut self.id().to_be_bytes().to_vec());

        // 4 bytes for Length of Disease Name
        res.append(&mut usize_to_u32(name_length).to_be_bytes().to_vec());

        // Disease name (n bytes)
        for c in name {
            res.push(*c);
        }

        // 4 bytes for number of HPO terms
        res.append(&mut usize_to_u32(self.hpo_terms().len()).to_be_bytes().to_vec());

        // HPO terms
        res.append(&mut self.hpo_terms().as_bytes());

        res
    }

    /// Returns an [`Disease`] from a bytes vector
    ///
    /// The byte layout for this method is defined in
    /// [`Disease::as_bytes`]
    ///
    /// # Errors
    ///
    /// This method can fail with a `ParseBinaryError` if the input
    /// bytes don't have the correct length or the name is not UTF-8
    /// encoded.
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
    fn from_bytes(bytes: &[u8]) -> Result<Self, HpoError> {
        // minimum length for a Disease without name and no HPO terms
        // This check is important because we're accessing the bytes
        // for size and ID directly and don't want to panic
        if bytes.len() < 4 + 4 + 4 + 4 {
            error!("Too few bytes for a Disease");
            return Err(HpoError::ParseBinaryError);
        }

        let total_len = u32_from_bytes(&bytes[0..]) as usize;

        if bytes.len() != total_len {
            error!(
                "Too few bytes to build a Disease. Expected {}, received {}",
                total_len,
                bytes.len()
            );
            return Err(HpoError::ParseBinaryError);
        }

        let id = u32_from_bytes(&bytes[4..]);
        let name_len = u32_from_bytes(&bytes[8..]) as usize;

        // Minimum length considering the name
        if bytes.len() < 16 + name_len {
            error!("Too few bytes for a Disease (including the name)");
            return Err(HpoError::ParseBinaryError);
        }

        let Ok(name) = String::from_utf8(bytes[12..12 + name_len].to_vec()) else {
            error!("Unable to parse the name of the Disease. It should consist of only utf8 characters");
            return Err(HpoError::ParseBinaryError);
        };

        let mut disease = Self::new(id.into(), &name);

        // An index for accessing the bytes during the iteration
        let mut idx_terms = 12 + name_len;
        let n_terms = u32_from_bytes(&bytes[idx_terms..]);

        if bytes.len() < 16 + name_len + n_terms as usize * 4 {
            error!(
                "Too few bytes in {}. {} terms, but {} bytes",
                name,
                n_terms,
                bytes.len()
            );
            return Err(HpoError::ParseBinaryError);
        }

        idx_terms += 4;
        for _ in 0..n_terms {
            let term_id = HpoTermId::from([
                bytes[idx_terms],
                bytes[idx_terms + 1],
                bytes[idx_terms + 2],
                bytes[idx_terms + 3],
            ]);
            idx_terms += 4;

            disease.add_term(term_id);
        }

        if idx_terms == total_len && idx_terms == bytes.len() {
            Ok(disease)
        } else {
            error!(
                "The length of the bytes blob did not match: {} vs {}",
                total_len, idx_terms
            );
            Err(HpoError::ParseBinaryError)
        }
    }
}

/// [`Disease`] Iterator
pub struct DiseaseIterator<'a, DID> {
    pub(crate) ontology: &'a Ontology,
    pub(crate) diseases: std::collections::hash_set::Iter<'a, DID>,
}

impl<'a, DID> DiseaseIterator<'a, DID> {
    /// Initialize a new [`DiseaseIterator`]
    ///
    /// This method requires the [`Ontology`] as a parameter since
    /// the actual [`Disease`] entities are stored in it.
    pub fn new(diseases: &'a HashSet<DID>, ontology: &'a Ontology) -> Self {
        DiseaseIterator {
            diseases: diseases.iter(),
            ontology,
        }
    }
}

impl<DID> Debug for DiseaseIterator<'_, DID> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "DiseaseIterator")
    }
}
