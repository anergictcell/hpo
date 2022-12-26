use crate::annotations::{GeneId, Genes};
use crate::annotations::{OmimDiseaseId, OmimDiseases};
use crate::term::HpoGroup;
use crate::term::InformationContent;
use crate::term::{HpoChildren, HpoParents, HpoTermId};
use crate::DEFAULT_NUM_OMIM;
use crate::DEFAULT_NUM_PARENTS;
use crate::{HpoError, DEFAULT_NUM_GENES};
use crate::{HpoResult, DEFAULT_NUM_ALL_PARENTS};

#[derive(Debug)]
pub(crate) struct HpoTermInternal {
    id: HpoTermId,
    name: String,
    parents: HpoParents,
    all_parents: HpoParents,
    children: HpoChildren,
    genes: Genes,
    omim_diseases: OmimDiseases,
    ic: InformationContent,
    #[allow(dead_code)]
    obsolete: bool,
    #[allow(dead_code)]
    replacement: Option<HpoTermId>,
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
            children: HpoChildren::with_capacity(DEFAULT_NUM_PARENTS),
            genes: Genes::with_capacity(DEFAULT_NUM_GENES),
            omim_diseases: OmimDiseases::with_capacity(DEFAULT_NUM_OMIM),
            ic: InformationContent::default(),
            obsolete: false,
            replacement: None,
        }
    }

    pub fn try_new(id: &str, name: &str) -> HpoResult<HpoTermInternal> {
        let id = HpoTermId::try_from(id)?;
        Ok(HpoTermInternal {
            id,
            name: name.to_string(),
            parents: HpoGroup::with_capacity(DEFAULT_NUM_PARENTS),
            all_parents: HpoGroup::with_capacity(DEFAULT_NUM_ALL_PARENTS),
            children: HpoChildren::with_capacity(DEFAULT_NUM_PARENTS),
            genes: Genes::with_capacity(DEFAULT_NUM_GENES),
            omim_diseases: OmimDiseases::with_capacity(DEFAULT_NUM_OMIM),
            ic: InformationContent::default(),
            obsolete: false,
            replacement: None,
        })
    }

    pub fn id(&self) -> &HpoTermId {
        &self.id
    }

    pub fn name(&self) -> &str {
        &self.name
    }

    pub fn parents(&self) -> &HpoParents {
        &self.parents
    }

    pub fn children(&self) -> &HpoChildren {
        &self.children
    }

    pub fn all_parents(&self) -> &HpoParents {
        &self.all_parents
    }

    pub fn all_parents_mut(&mut self) -> &mut HpoParents {
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

    pub fn add_parent(&mut self, parent_id: HpoTermId) {
        self.parents.insert(parent_id);
    }

    pub fn add_child(&mut self, child_id: HpoTermId) {
        self.children.insert(child_id);
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

    #[allow(dead_code)]
    pub fn obsolete(&self) -> bool {
        self.obsolete
    }

    #[allow(dead_code)]
    pub fn obsolete_mut(&mut self) -> &mut bool {
        &mut self.obsolete
    }

    #[allow(dead_code)]
    pub fn replacement(&self) -> Option<HpoTermId> {
        self.replacement
    }

    #[allow(dead_code)]
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
    ///
    pub fn as_bytes(&self) -> Vec<u8> {
        // 4 bytes for total length
        // 4 bytes for TermID (big-endian)
        // 1 byte for Name length (u8) -> Name cannot be longer than 255 bytes
        // name in u8 encoded
        let name = self.name().as_bytes();
        let name_length = std::cmp::min(name.len(), 255);
        let size = name_length + 4 + 4 + 1;

        let mut res = Vec::with_capacity(size);

        // 4 bytes for total length
        res.append(&mut (size as u32).to_be_bytes().to_vec());

        // 4 bytes to Term-ID
        res.append(&mut self.id.to_be_bytes().to_vec());

        // 1 byte for Length of Term Name (can't be longer than 255 bytes)
        res.push(name_length as u8);

        // Term name (up to 255 bytes)
        for c in name.iter().take(255) {
            res.push(*c);
        }
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
    pub fn parents_as_byte(&self) -> Vec<u8> {
        let mut term_parents: Vec<u8> = Vec::new();
        let n_parents = self.parents().len() as u32;
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

impl TryFrom<&[u8]> for HpoTermInternal {
    type Error = HpoError;
    /// Crates an `HpoTermInternal` from raw bytes
    ///
    /// See [`HpoTermInternal::as_bytes`] for description of the byte layout
    fn try_from(bytes: &[u8]) -> Result<Self, Self::Error> {
        if bytes.len() < 4 + 4 + 1 {
            return Err(HpoError::ParseBinaryError);
        }
        let total_len = u32::from_be_bytes([bytes[0], bytes[1], bytes[2], bytes[3]]);

        let id = u32::from_be_bytes([bytes[4], bytes[5], bytes[6], bytes[7]]);
        let name_len = bytes[8] as usize;

        if bytes.len() < 4 + 4 + 1 + name_len {
            return Err(HpoError::ParseBinaryError);
        }

        let name = match String::from_utf8(bytes[9..total_len as usize].to_vec()) {
            Ok(s) => s,
            Err(_) => return Err(HpoError::ParseBinaryError),
        };
        Ok(HpoTermInternal::new(name, id.into()))
    }
}

/// Builder to crate multiple [`HpoTermInternal`] from raw bytes
pub(crate) struct BinaryTermBuilder<'a> {
    bytes: &'a [u8],
    idx: usize,
}

impl<'a> BinaryTermBuilder<'a> {
    /// Crates a new [`BinaryTermBuilder`]
    pub fn new(bytes: &'a [u8]) -> Self {
        BinaryTermBuilder { bytes, idx: 0 }
    }
}

impl Iterator for BinaryTermBuilder<'_> {
    type Item = HpoTermInternal;
    fn next(&mut self) -> Option<Self::Item> {
        let bytes = &self.bytes[self.idx..];

        if bytes.is_empty() {
            return None;
        }

        let term_len = u32::from_be_bytes(bytes[0..4].try_into().unwrap()) as usize;

        assert!(
            (bytes.len() >= term_len),
            "Invalid bytes left over in BinaryTermBuilder"
        );

        self.idx += term_len;
        let term = &bytes[..term_len];
        Some(HpoTermInternal::try_from(term).unwrap())
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn to_bytes() {
        let term = HpoTermInternal::new(String::from("Foobar"), 123u32.into());

        let bytes: Vec<u8> = term.as_bytes();

        let term_len = u32::from_be_bytes(bytes[0..4].try_into().unwrap()) as usize;
        assert_eq!(term_len, 4 + 4 + 1 + 6);
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
        let term2 = HpoTermInternal::try_from(&bytes[..]).unwrap();
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
            v.append(&mut t.as_bytes());
        }

        let mut term_iter = BinaryTermBuilder::new(&v);

        for (name, id) in test_terms {
            let term = term_iter.next().unwrap();
            assert_eq!(term.name(), name);
            assert_eq!(term.id().as_u32(), id);
        }

        assert!(term_iter.next().is_none());
    }
}
