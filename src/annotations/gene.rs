use core::fmt::Debug;
use std::cmp::PartialEq;
use std::collections::HashSet;
use std::convert::TryFrom;
use std::fmt::Display;
use std::hash::Hash;

use tracing::error;

use crate::annotations::AnnotationId;
use crate::set::HpoSet;
use crate::term::HpoGroup;
use crate::u32_from_bytes;
use crate::HpoError;
use crate::HpoResult;
use crate::HpoTermId;
use crate::Ontology;

/// A set of genes
///
/// The set does not contain [`Gene`]s itself, but only their [`GeneId`]s.
/// Currently implemented using [`HashSet`] but any other implementation
/// should work as well given that each [`GeneId`] must appear only once
/// and it provides an iterator of [`GeneId`]
pub type Genes = HashSet<GeneId>;

/// A unique identifier for a [`Gene`]
///
/// This value can - in theory - represent any numerical unique value.
/// When using the default JAX provided masterdata, it represents
/// the NCBI Gene ID
#[derive(Clone, Copy, Default, Debug, Hash, PartialEq, PartialOrd, Eq, Ord)]
pub struct GeneId {
    inner: u32,
}

impl AnnotationId for GeneId {
    /// Convert `self` to `u32`
    fn as_u32(&self) -> u32 {
        self.inner
    }
}

impl TryFrom<&str> for GeneId {
    type Error = HpoError;
    fn try_from(value: &str) -> Result<Self, Self::Error> {
        Ok(GeneId {
            inner: value.parse::<u32>()?,
        })
    }
}

impl From<u32> for GeneId {
    fn from(inner: u32) -> Self {
        GeneId { inner }
    }
}

impl Display for GeneId {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "NCBI-GeneID:{}", self.inner)
    }
}

/// A single gene
///
/// A gene has a unique [`GeneId`] and a name (symbol) and is
/// connected to a set of HPO terms
#[derive(Default, Debug, Clone)]
pub struct Gene {
    id: GeneId,
    name: String,
    hpos: HpoGroup,
}

impl Gene {
    /// Initializes a new Gene
    ///
    /// This method should rarely, if ever, be used directly. The
    /// preferred way to create new genes is through [`Ontology::add_gene`]
    /// to ensure that each gene exists only once.
    pub fn new(id: GeneId, name: &str) -> Gene {
        Gene {
            id,
            name: name.to_string(),
            hpos: HpoGroup::default(),
        }
    }

    /// Initializes a new Gene from `str` values
    ///
    /// This method should rarely, if ever, be used directly. The
    /// preferred way to create new genes is through [`Ontology::add_gene`]
    /// to ensure that each gene exists only once.
    ///
    /// # Errors
    ///
    /// If the id is not a correct `GeneId`, returns [`HpoError::ParseIntError`]
    pub fn from_parts(id: &str, name: &str) -> HpoResult<Gene> {
        Ok(Gene {
            id: GeneId::try_from(id)?,
            name: name.to_string(),
            hpos: HpoGroup::default(),
        })
    }

    /// The unique [`GeneId`] of the gene, most likely the NCBI Gene ID
    pub fn id(&self) -> &GeneId {
        &self.id
    }

    /// The name of the gene (gene symbol)
    pub fn name(&self) -> &str {
        &self.name
    }

    /// The gene symbol (identical to [`Gene::id`])
    pub fn symbol(&self) -> &str {
        &self.name
    }

    /// The set of connected HPO terms
    pub fn hpo_terms(&self) -> &HpoGroup {
        &self.hpos
    }

    /// Returns a binary representation of the `Gene`
    ///
    /// The binary layout is defined as:
    ///
    /// | Byte offset | Number of bytes | Description |
    /// | --- | --- | --- |
    /// | 0 | 4 | The total length of the binary data blob as big-endian `u32` |
    /// | 4 | 4 | The Gene ID as big-endian `u32` |
    /// | 8 | 1 | The length of the Gene Name / Symbol (converted to a u8 vector) as a `u8` |
    /// | 9 | n | The Gene name/symbol as u8 vector. If the name has more than 255 bytes, it is trimmed to 255 |
    /// | 9 + n | 4 | The number of associated HPO terms as big-endian `u32` |
    /// | 13 + n | x * 4 | The HPO Term IDs of the associated terms, each encoded as big-endian `u32` |
    ///
    /// # Examples
    ///
    /// ```
    /// use hpo::annotations::Gene;
    ///
    /// let mut gene = Gene::from_parts("123", "FooBar").unwrap();
    /// let bytes = gene.as_bytes();
    ///
    /// assert_eq!(bytes.len(), 4 + 4 + 1 + 6 + 4);
    /// assert_eq!(bytes[4..8], [0u8, 0u8, 0u8, 123u8]);
    /// assert_eq!(bytes[8], 6u8);
    /// ```
    pub fn as_bytes(&self) -> Vec<u8> {
        fn usize_to_u32(n: usize) -> u32 {
            n.try_into().expect("unable to convert {n} to u32")
        }
        let name = self.name().as_bytes();
        let name_length = std::cmp::min(name.len(), 255);
        let size = 4 + 4 + 1 + name_length + 4 + self.hpos.len() * 4;

        let mut res = Vec::new();

        // 4 bytes for total length
        res.append(&mut usize_to_u32(size).to_be_bytes().to_vec());

        // 4 bytes for Gene-ID
        res.append(&mut self.id.to_be_bytes().to_vec());

        // 1 byte for Length of Gene Name (can't be longer than 255 bytes)
        // casting is safe, since name_length is < 256
        #[allow(clippy::cast_possible_truncation)]
        res.push(name_length as u8);

        // Gene name/symbol (up to 255 bytes)
        for c in name.iter().take(255) {
            res.push(*c);
        }

        // 4 bytes for number of HPO terms
        res.append(&mut usize_to_u32(self.hpos.len()).to_be_bytes().to_vec());

        // HPO terms
        res.append(&mut self.hpos.as_bytes());

        res
    }

    /// Returns an [`HpoSet`] from the `Gene`
    pub fn to_hpo_set<'a>(&self, ontology: &'a Ontology) -> HpoSet<'a> {
        HpoSet::new(ontology, self.hpos.clone())
    }

    /// Connect another [HPO term](`crate::HpoTerm`) to the gene
    ///
    /// # Note
    ///
    /// This method does **not** add the [`Gene`] to the [HPO term](`crate::HpoTerm`).
    /// Clients should not use this method, unless they are creating their own Ontology.
    pub fn add_term<I: Into<HpoTermId>>(&mut self, term_id: I) -> bool {
        self.hpos.insert(term_id)
    }
}

impl PartialEq for Gene {
    fn eq(&self, other: &Gene) -> bool {
        self.id == other.id
    }
}
impl Eq for Gene {}

impl TryFrom<&[u8]> for Gene {
    type Error = HpoError;
    /// Returns a [`Gene`] from a bytes vector
    ///
    /// The byte layout for this method is defined in
    /// [`Gene::as_bytes`]
    ///
    /// # Examples
    ///
    /// ```
    /// use hpo::annotations::{Gene, GeneId};
    ///
    /// let bytes = vec![
    ///     0u8, 0u8, 0u8, 19u8, // Total size of Blop
    ///     0u8, 0u8, 0u8, 123u8, // ID of the gene => 123
    ///     6u8, // Length of gene symbol
    ///     b'F', b'o', b'o', b'b', b'a', b'r', // Foobar
    ///     0u8, 0u8, 0u8, 0u8  // Number of associated HPO Terms => 0
    /// ];
    /// let gene = Gene::try_from(&bytes[..]).unwrap();
    ///
    /// assert_eq!(gene.name(), "Foobar");
    /// assert_eq!(gene.id(), &GeneId::from(123u32));
    /// ```
    fn try_from(bytes: &[u8]) -> Result<Self, Self::Error> {
        // minimum length for a Gene without name and no HPO terms
        // This check is important because we're accessing the bytes
        // for size and ID directly and don't want to panic
        if bytes.len() < 4 + 4 + 1 + 4 {
            error!("Too few bytes for a Gene");
            return Err(HpoError::ParseBinaryError);
        }
        let total_len = u32_from_bytes(&bytes[0..]) as usize;

        if bytes.len() != total_len {
            error!(
                "Too few bytes to build gene. Expected {}, received {}",
                total_len,
                bytes.len()
            );
            return Err(HpoError::ParseBinaryError);
        }

        let id = u32_from_bytes(&bytes[4..]);
        let name_len = bytes[8] as usize;

        // Minimum length considering the name
        if bytes.len() < 13 + name_len {
            error!("Too few bytes for an Gene (including the name)");
            return Err(HpoError::ParseBinaryError);
        }

        let Ok(name) = String::from_utf8(bytes[9..9 + name_len].to_vec()) else {
            error!("Unable to parse the name of the Gene");
            return Err(HpoError::ParseBinaryError);
        };

        let mut gene = Gene::new(id.into(), &name);

        let mut idx_terms = 9 + name_len;
        let n_terms = u32_from_bytes(&bytes[idx_terms..]);

        if bytes.len() < 13 + name_len + n_terms as usize * 4 {
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

            gene.add_term(term_id);
        }

        if idx_terms == total_len && idx_terms == bytes.len() {
            Ok(gene)
        } else {
            error!(
                "The length of the bytes blob did not match: {} vs {}",
                total_len, idx_terms
            );
            Err(HpoError::ParseBinaryError)
        }
    }
}

impl Hash for Gene {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.id.hash(state);
    }
}

/// [`Gene`] Iterator
pub struct GeneIterator<'a> {
    ontology: &'a Ontology,
    genes: std::collections::hash_set::Iter<'a, GeneId>,
}

impl<'a> GeneIterator<'a> {
    /// Initialize a new [`GeneIterator`]
    ///
    /// This method requires the [`Ontology`] as a parameter since
    /// the actual [`Gene`] entities are stored in it and not in [`Genes`]
    /// itself
    pub fn new(genes: &'a Genes, ontology: &'a Ontology) -> Self {
        GeneIterator {
            genes: genes.iter(),
            ontology,
        }
    }
}

impl<'a> std::iter::Iterator for GeneIterator<'a> {
    type Item = &'a Gene;
    fn next(&mut self) -> Option<Self::Item> {
        match self.genes.next() {
            Some(gene_id) => Some(self.ontology.gene(gene_id).unwrap()),
            None => None,
        }
    }
}

impl Debug for GeneIterator<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "GeneIterator")
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn gene_to_binary() {
        let mut gene = Gene::from_parts("123", "FooBar").unwrap();
        gene.add_term(66u32);
        gene.add_term(77u32);

        let bin = gene.as_bytes();

        assert_eq!(bin.len(), 4 + 4 + 1 + 6 + 4 + 8);
    }

    #[test]
    fn gene_to_and_from_binary() {
        let mut gene = Gene::from_parts("123", "FooBar").unwrap();
        gene.add_term(66u32);
        gene.add_term(77u32);

        let bin = gene.as_bytes();

        let gene2: Gene = bin[..].try_into().unwrap();

        assert_eq!(gene.name(), gene2.name());
        assert_eq!(gene.id(), gene2.id());
        assert_eq!(gene.hpo_terms().len(), gene2.hpo_terms().len());
        for (a, b) in gene.hpo_terms().iter().zip(gene2.hpo_terms().iter()) {
            assert_eq!(a, b);
        }
    }

    #[test]
    fn gene_with_non_utf8_name() {
        let mut gene = Gene::new(123u32.into(), "Foo😀Bar");
        gene.add_term(66u32);
        gene.add_term(77u32);

        let bin = gene.as_bytes();

        // the smiley uses 4 bytes, so the name length is 10
        assert_eq!(bin.len(), 4 + 4 + 1 + 10 + 4 + 8);

        let gene2 = Gene::try_from(&bin[..]).expect("Can't build Gene");

        assert_eq!(gene.name(), gene2.name());
        assert_eq!(gene.id(), gene2.id());
        assert_eq!(gene.hpo_terms().len(), gene2.hpo_terms().len());
        for (a, b) in gene.hpo_terms().iter().zip(gene2.hpo_terms().iter()) {
            assert_eq!(a, b);
        }
    }

    #[test]
    fn gene_with_long_name() {
        let name = "Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat. Duis aute irure dolor in reprehenderit in voluptate velit esse cillum dolore eu fugiat nulla pariatur. Excepteur sint occaecat cupidatat non proident, sunt in culpa qui officia deserunt mollit anim id est laborum.";
        let mut gene = Gene::new(123u32.into(), name);
        gene.add_term(66u32);
        gene.add_term(77u32);

        let bin = gene.as_bytes();

        let gene2 = Gene::try_from(&bin[..]).expect("Can't build Disease");

        assert_ne!(gene.name(), gene2.name());
        assert_eq!(gene.name()[0..255], gene2.name()[0..255]);
        assert_eq!(gene.id(), gene2.id());
        assert_eq!(gene.hpo_terms().len(), gene2.hpo_terms().len());
        for (a, b) in gene.hpo_terms().iter().zip(gene2.hpo_terms().iter()) {
            assert_eq!(a, b);
        }
    }

    #[test]
    fn gene_with_wrong_length() {
        let mut gene = Gene::new(123u32.into(), "foobar");
        gene.add_term(66u32);
        gene.add_term(77u32);

        let mut bin = gene.as_bytes();

        assert!(Gene::try_from(&bin[..12]).is_err());
        assert!(Gene::try_from(&bin[..18]).is_err());
        assert!(Gene::try_from(&bin[..26]).is_err());
        assert!(Gene::try_from(&bin[..27]).is_ok());

        bin.push(1);
        assert!(Gene::try_from(&bin[..27]).is_ok());
        assert!(Gene::try_from(&bin[..28]).is_err());
    }
}
