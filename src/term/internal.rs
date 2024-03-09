use crate::annotations::AnnotationId;
use crate::parser::binary::term::{from_bytes_v1, from_bytes_v2};
use crate::parser::binary::{BinaryVersion, Bytes};
use std::hash::Hash;

use crate::annotations::{GeneId, Genes};
use crate::annotations::{OmimDiseaseId, OmimDiseases};
use crate::term::{HpoGroup, HpoTermId, InformationContent};
use crate::DEFAULT_NUM_PARENTS;
use crate::{HpoError, DEFAULT_NUM_GENES};
use crate::{HpoResult, DEFAULT_NUM_ALL_PARENTS};
use crate::{HpoTerm, DEFAULT_NUM_OMIM};

#[derive(Clone, Debug)]
pub(crate) struct HpoTermInternal {
    id: HpoTermId,
    name: String,
    parents: HpoGroup,
    all_parents: HpoGroup,
    children: HpoGroup,
    genes: Genes,
    omim_diseases: OmimDiseases,
    ic: InformationContent,
    obsolete: bool,
    replacement: Option<HpoTermId>,
}

impl Hash for HpoTermInternal {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.id.hash(state);
    }
}

impl Default for HpoTermInternal {
    fn default() -> Self {
        HpoTermInternal::new(String::from("HP:0000000"), 0u32.into())
    }
}

impl HpoTermInternal {
    pub fn new(name: String, id: HpoTermId) -> HpoTermInternal {
        HpoTermInternal {
            id,
            name,
            parents: HpoGroup::with_capacity(DEFAULT_NUM_PARENTS),
            all_parents: HpoGroup::with_capacity(DEFAULT_NUM_ALL_PARENTS),
            children: HpoGroup::with_capacity(DEFAULT_NUM_PARENTS),
            genes: Genes::with_capacity(DEFAULT_NUM_GENES),
            omim_diseases: OmimDiseases::with_capacity(DEFAULT_NUM_OMIM),
            ic: InformationContent::default(),
            obsolete: false,
            replacement: None,
        }
    }

    pub fn try_new(id: &str, name: &str) -> HpoResult<HpoTermInternal> {
        let id = HpoTermId::try_from(id)?;
        Ok(Self::new(name.to_string(), id))
    }

    pub fn id(&self) -> &HpoTermId {
        &self.id
    }

    pub fn name(&self) -> &str {
        &self.name
    }

    pub fn parents(&self) -> &HpoGroup {
        &self.parents
    }

    pub fn children(&self) -> &HpoGroup {
        &self.children
    }

    pub fn all_parents(&self) -> &HpoGroup {
        &self.all_parents
    }

    pub fn all_parents_mut(&mut self) -> &mut HpoGroup {
        &mut self.all_parents
    }

    pub fn genes(&self) -> &Genes {
        &self.genes
    }

    pub fn omim_diseases(&self) -> &OmimDiseases {
        &self.omim_diseases
    }

    pub fn parents_cached(&self) -> bool {
        if self.parents.is_empty() {
            true
        } else {
            !self.all_parents.is_empty()
        }
    }

    pub fn add_parent<I: Into<HpoTermId>>(&mut self, parent_id: I) {
        self.parents.insert(parent_id.into());
    }

    pub fn add_child<I: Into<HpoTermId>>(&mut self, child_id: I) {
        self.children.insert(child_id.into());
    }

    pub fn add_gene(&mut self, gene_id: GeneId) -> bool {
        self.genes.insert(gene_id)
    }

    pub fn add_omim_disease(&mut self, omim_disease_id: OmimDiseaseId) -> bool {
        self.omim_diseases.insert(omim_disease_id)
    }

    pub fn information_content(&self) -> &InformationContent {
        &self.ic
    }

    pub fn information_content_mut(&mut self) -> &mut InformationContent {
        &mut self.ic
    }

    pub fn obsolete(&self) -> bool {
        self.obsolete
    }

    pub fn obsolete_mut(&mut self) -> &mut bool {
        &mut self.obsolete
    }

    pub fn replacement(&self) -> Option<HpoTermId> {
        self.replacement
    }

    pub fn replacement_mut(&mut self) -> &mut Option<HpoTermId> {
        &mut self.replacement
    }

    /// Returns a binary representation of the `HpoTermInternal`
    ///
    /// The binary layout is defined as:
    ///
    /// | Byte offset | Number of bytes | Description |
    /// | --- | --- | --- |
    /// | 0 | 4 | The total length of the binary data blob as big-endian `u32` |
    /// | 4 | 4 | The Term ID as big-endian `u32` |
    /// | 8 | 1 | The length of the Term Name (converted to a u8 vector) as a `u8` |
    /// | 9 | n | The Term name as u8 vector. If the name has more than 255 bytes, it is trimmed to 255 |
    /// | 9 + n | 1 | Flag to indicate if term is obsolete
    /// | 10 + n | 4 | Term ID of a replacement term as big-endian `u32` or `0` if `None` |
    ///
    /// # Panics
    ///
    /// This method will panic if the total byte length is longer than `u32::MAX`
    pub fn as_bytes(&self) -> Vec<u8> {
        // 4 bytes for total length
        // 4 bytes for TermID (big-endian)
        // 1 byte for Name length (u8) -> Name cannot be longer than 255 bytes
        // 1 byte for obsolete flag
        // 4 byte for replacement term
        // name in u8 encoded
        let name = self.name().as_bytes();
        let name_length = std::cmp::min(name.len(), 255);
        let size = name_length + 4 + 4 + 1 + 1 + 4;

        let mut res = Vec::with_capacity(size);

        // 4 bytes for total length
        res.append(&mut u32::try_from(size).unwrap().to_be_bytes().to_vec());

        // 4 bytes to Term-ID
        res.append(&mut self.id.to_be_bytes().to_vec());

        // 1 byte for Length of Term Name (can't be longer than 255 bytes)
        // casting is safe, since name_length is < 256
        #[allow(clippy::cast_possible_truncation)]
        res.push(name_length as u8);

        // Term name (up to 255 bytes)
        for c in name.iter().take(name_length) {
            res.push(*c);
        }

        // 1 byte for various flags, currently only obsolete flag
        if self.obsolete {
            res.push(1u8);
        } else {
            res.push(0u8);
        }

        // 4 bytes for replace term (or 0 if `None`)
        res.append(
            &mut self
                .replacement
                .unwrap_or(0u32.into())
                .to_be_bytes()
                .to_vec(),
        );

        res
    }

    /// Returns a binary representation of Term - Parent connections
    ///
    /// The binary layout is defined as:
    ///
    /// | Byte offset | Number of bytes | Description |
    /// | --- | --- | --- |
    /// | 0 | 4 | The number of parent terms as big-endian `u32` |
    /// | 4 | 4 | The Term ID of the term as big-endian `u32` |
    /// | 8 | 4 * n | The Term ID of all parents as big-endian `u32` |
    ///
    /// # Panics
    ///
    /// This method will panic if there are more than `u32::MAX` parents
    pub fn parents_as_byte(&self) -> Vec<u8> {
        let mut term_parents: Vec<u8> = Vec::new();
        let n_parents: u32 = self.parents().len().try_into().unwrap();
        term_parents.append(&mut n_parents.to_be_bytes().to_vec());
        term_parents.append(&mut self.id().to_be_bytes().to_vec());
        for parent in self.parents() {
            term_parents.append(&mut parent.to_be_bytes().to_vec());
        }
        term_parents
    }
}

impl PartialEq for HpoTermInternal {
    fn eq(&self, other: &Self) -> bool {
        self.id == other.id
    }
}

impl Eq for HpoTermInternal {}

impl TryFrom<Bytes<'_>> for HpoTermInternal {
    type Error = HpoError;
    /// Crates an `HpoTermInternal` from raw bytes
    ///
    /// See [`HpoTermInternal::as_bytes`] for description of the byte layout
    fn try_from(bytes: Bytes) -> Result<Self, Self::Error> {
        match bytes.version() {
            BinaryVersion::V1 => from_bytes_v1(bytes),
            BinaryVersion::V2 => from_bytes_v2(bytes),
        }
    }
}

impl<'a> From<&HpoTerm<'a>> for HpoTermInternal {
    fn from(term: &HpoTerm) -> Self {
        let mut internal = Self::new(term.name().to_string(), term.id());
        *internal.obsolete_mut() = term.is_obsolete();
        *internal.replacement_mut() = term.replaced_by().map(|repl| repl.id());
        internal
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::parser::binary::BinaryTermBuilder;

    #[test]
    fn to_bytes() {
        let term = HpoTermInternal::new(String::from("Foobar"), 123u32.into());

        let bytes: Vec<u8> = term.as_bytes();

        let term_len = u32::from_be_bytes(bytes[0..4].try_into().unwrap()) as usize;
        assert_eq!(term_len, 4 + 4 + 1 + 6 + 5);
        let term_id = u32::from_be_bytes(bytes[4..8].try_into().unwrap());
        assert_eq!(term_id, 123);
        let name_len = bytes[8] as usize;
        assert_eq!(name_len, 6);
        let name = String::from_utf8(bytes[9..9 + name_len].to_vec()).unwrap();
        assert_eq!(name, "Foobar");
    }

    #[test]
    fn from_bytes() {
        let term = HpoTermInternal::new(String::from("Foobar"), 123u32.into());
        let bytes: Vec<u8> = term.as_bytes();
        let term2 = HpoTermInternal::try_from(Bytes::new(&bytes[..], BinaryVersion::V2)).unwrap();
        assert_eq!(term2.name(), term.name());
        assert_eq!(term2.id(), term.id());
    }

    #[test]
    fn from_multiple_bytes() {
        let mut v: Vec<u8> = Vec::new();

        let test_terms = [
            ("t1", 1u32),
            ("Term with a very long name", 2u32),
            ("", 3u32),
            ("Abnormality", 4u32),
        ];

        for (name, id) in test_terms {
            let t = HpoTermInternal::new(String::from(name), id.into());
            println!("Building: {t:?}");
            v.append(&mut t.as_bytes());
        }

        let mut term_iter = BinaryTermBuilder::new(Bytes::new(&v, BinaryVersion::V2));

        for (name, id) in test_terms {
            let term = term_iter.next().unwrap();
            println!("Checking: {term:?} [{name}-{id}]");
            assert_eq!(term.name(), name);
            assert_eq!(term.id().as_u32(), id);
        }

        assert!(term_iter.next().is_none());
    }
}
