use core::fmt::Debug;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::Read;
use std::ops::BitOr;
use std::path::Path;

use tracing::debug;

use crate::annotations::AnnotationId;
use crate::annotations::{Gene, GeneId};
use crate::annotations::{OmimDisease, OmimDiseaseId};
use crate::parser;
use crate::parser::binary::{BinaryTermBuilder, BinaryVersion, Bytes};
use crate::term::internal::HpoTermInternal;
use crate::term::{HpoGroup, HpoTerm};
use crate::u32_from_bytes;
use crate::HpoResult;
use crate::{HpoError, HpoTermId};

pub mod comparison;
mod termarena;
use comparison::Comparison;
use termarena::Arena;

#[cfg_attr(doc, aquamarine::aquamarine)]
/// `Ontology` is the main interface of the `hpo` crate and contains all data
///
/// The [`Ontology`] struct holds all information about the ontology
/// and all [`HpoTerm`]s, [`Gene`]s and [`OmimDisease`]s.
///
/// # Examples
///
/// ```
/// use hpo::{Ontology, HpoTermId};
///
/// let ontology = Ontology::from_binary("tests/example.hpo").unwrap();
///
/// // get single terms from the ontology
///
/// let absent_term = HpoTermId::try_from("HP:9999999").unwrap();
/// assert!(ontology.hpo(absent_term).is_none());
///
/// let present_term = HpoTermId::try_from("HP:0000001").unwrap();
/// let root_term = ontology.hpo(present_term).unwrap();
/// assert_eq!(root_term.name(), "All");
///
/// // simplified way to get an `HpoTerm` by using the `u32` `HpoTermId`
/// let term = ontology.hpo(118u32).unwrap();
/// assert_eq!(term.name(), "Phenotypic abnormality");
///
/// // get all genes of the ontology
/// assert_eq!(ontology.genes().count(), 4852);
///
/// // get all diseases of the ontology
/// assert_eq!(ontology.omim_diseases().count(), 4431);
///
/// // Iterate all HPO terms
/// for term in &ontology {
///     // do something with term
///     println!("{}", term.name());
/// }
/// ```
///
/// # Construction
///
/// There are two main ways to build the Ontology
/// 1. Download the standard annotation data from
///     [Jax HPO](https://hpo.jax.org/) itself.
///     Then use [`Ontology::from_standard`] to load the data.
///     You need the following files:
///     - `phenotype.hpoa` (Required to connect [`OmimDisease`]s to [`HpoTerm`]s)
///     - `phenotype_to_genes.txt` (Required to connect [`Gene`]s to [`HpoTerm`]s)
///     - `hp.obo` (Required for [`HpoTerm`]s and their connection to each other)
/// 2. Load the ontology from a binary build using [`Ontology::from_binary`].
///
///     The [Github repository](https://github.com/anergictcell/hpo) of this crate
///     contains a binary build of the ontology
///     <https://github.com/anergictcell/hpo/blob/main/tests/ontology.hpo>.
///     The snapshot will not always be up to date, so please double-check yourself.
///
///     You can crate your own binary build of the ontology using the
///     `examples/obo_to_bin.rs` example.
///
///     `cargo run --example --release obo_to_bin <PATH TO FOLDER WITH JAX DATA> <OUTPUT FILENAME>`
///
/// You can also build it all by yourself (not recommended), in which case you
/// will have to:
/// 1. construct an empty Ontology [`Ontology::default`]
/// 2. Add all terms [`Ontology::insert_term`]
/// 3. Connect terms to their parents [`Ontology::add_parent`]
/// 4. Cache all parent, child and grandparent connections [`Ontology::create_cache`]
/// 5. Add genes and diseases to the ontology
///     - [`Ontology::add_gene`] and [`Ontology::add_omim_disease`]
///     - Connect genes and diseases to the [`HpoTerm`]s using
///         [`Ontology::link_gene_term`] and [`Ontology::link_omim_disease_term`]
///         (this will automatically take care of "inheriting" the connection to all
///         parent terms)
///     - make sure to also add the linked terms to the genes and diseases
///         [`Gene::add_term`] and [`OmimDisease::add_term`]
/// 6. Calculate the information content [`Ontology::calculate_information_content`]
///
///
/// # Layout
///
/// The [`Ontology`] contains all terms and all associated genes and diseases.
/// [`HpoTerm`]s are connected to each other in a directed relationship. Every term
/// (except the term `All`) has at least one parent term in an `is_a` relationship.
/// Terms and [`crate::annotations`] ([`Gene`]s, [`OmimDisease`]s) have a many-to-many relationship. The
/// [`Ontology`] does not contain a direct relationship between genes and diseases. This relation
/// is only present indirectly via the connected [`HpoTerm`]s.
///
/// ```mermaid
/// erDiagram
///     ONTOLOGY ||--|{ HPOTERM : contains
///     HPOTERM ||--|{ HPOTERM : is_a
///     HPOTERM }|--o{ DISEASE : phenotype_of
///     HPOTERM }|--o{ GENE : phenotype_of
///     HPOTERM {
///         str name
///         HpoTermId id
///         HpoTerms parents
///         HpoTerms children
///         Genes genes
///         OmimDiseases omim_diseases
///     }
///     DISEASE {
///         str name
///         OmimDiseaseId id
///         HpoGroup hpo_terms
///     }
///     GENE {
///         str name
///         GeneId id
///         HpoGroup hpo_terms
///     }
/// ```
///
/// # Relations of different public struct in this module
///
/// The below diagram looks complicated at first, but the
/// relationship of all entities follows a logical pattern.
/// `HpoTerm` and `HpoSet` are the most important public structs.
/// The `HpoGroup` is more relevant for internal use, but can also be
/// useful for fast set-based operations.
///
/// ```mermaid
/// classDiagram
///     class Ontology {
///         into_iter()
///     }
///
///     class HpoTerm{
///         - HpoTermId id
///         - &Ontology
///         parents() HpoTerms
///         parent_ids() HpoGroup
///         all_parent_ids() HpoGroup
///         children() HpoTerms
///         children_ids() HpoTerms
///         common_ancestors() Combine
///         union_ancestors() Combine
///         many-more()
///     }
///
///     class HpoGroup {
///         - Set~HpoTermId~
///         into_iter()
///         terms()
///     }
///
///     class HpoSet {
///         - HpoGroup
///         - &Ontology
///         similarity(...) f32
///         information_content()
///     }
///
///     class HpoTermId {
///         - u32: id
///     }
///
///     class `ontology::Iter` {
///         next() HpoTerm
///     }
///
///     class `term::Iter` {
///         next() HpoTerm
///     }
///
///     class `group::Iter` {
///         next() HpoTermId
///     }
///
///     class Combine {
///         - HpoGroup
///         into_iter()
///     }
///
///     Ontology ..|> `ontology::Iter`: hpos()
///     HpoSet ..|> `term::Iter`: iter()
///     HpoGroup ..|> `group::Iter`: iter()
///     HpoGroup ..|> `term::Iter`: terms()
///     Combine ..|> `term::Iter`: iter()
///
///     `ontology::Iter` --o HpoGroup: collect()
///     `ontology::Iter` --* HpoTerm: iterates()
///
///     `term::Iter` --* HpoTerm: iterates()
///     `term::Iter` --o HpoGroup: collect()
///
///     `group::Iter` --* HpoTermId: iterates()
///     `group::Iter` --o HpoGroup: collect()
///
///     HpoTerm ..|> HpoGroup: parent_ids()/children_ids()
///     HpoTerm ..|> `term::Iter`: parents()/children()
///     HpoTerm ..|> `Combine`: ..._ancestors()
/// ```
///
/// # Example ontology
///
/// For all examples and tests in this documentation, we're using the
/// following small subset of the full Ontology:
///
/// ```mermaid
/// graph TD
/// HP:0011017["HP:0011017<br>
/// Abnormal cellular physiology"]
/// HP:0010662["HP:0010662<br>
/// Abnormality of the diencephalon"]
/// HP:0010662 --> HP:0012285
/// HP:0000005["HP:0000005<br>
/// Mode of inheritance"]
/// HP:0000005 --> HP:0034345
/// HP:0012648["HP:0012648<br>
/// Decreased inflammatory response"]
/// HP:0012443["HP:0012443<br>
/// Abnormality of brain morphology"]
/// HP:0012443 --> HP:0100547
/// HP:0003674["HP:0003674<br>
/// Onset"]
/// HP:0003674 --> HP:0003581
/// HP:0010978["HP:0010978<br>
/// Abnormality of immune system physiology"]
/// HP:0010978 --> HP:0012647
/// HP:0000707["HP:0000707<br>
/// Abnormality of the nervous system"]
/// HP:0000707 --> HP:0012638
/// HP:0000707 --> HP:0012639
/// HP:0034345["HP:0034345<br>
/// Mendelian inheritance"]
/// HP:0034345 --> HP:0000007
/// HP:0000001["HP:0000001<br>
/// All"]
/// HP:0000001 -----> HP:0000005
/// HP:0000001 --> HP:0000118
/// HP:0000001 --> HP:0012823
/// HP:0000818["HP:0000818<br>
/// Abnormality of the endocrine system"]
/// HP:0000818 --> HP:0000864
/// HP:0100547["HP:0100547<br>
/// Abnormal forebrain morphology"]
/// HP:0100547 ----> HP:0010662
/// HP:0012647["HP:0012647<br>
/// Abnormal inflammatory response"]
/// HP:0012647 --> HP:0012648
/// HP:0001939["HP:0001939<br>
/// Abnormality of metabolism/homeostasis"]
/// HP:0001939 --> HP:0011017
/// HP:0001939 ---> HP:0025454
/// HP:0003581["HP:0003581<br>
/// Adult onset"]
/// HP:0012823["HP:0012823<br>
/// Clinical modifier"]
/// HP:0012823 --> HP:0031797
/// HP:0012285["HP:0012285<br>
/// Abnormal hypothalamus physiology"]
/// HP:0012638["HP:0012638<br>
/// Abnormal nervous system physiology"]
/// HP:0012638 ----> HP:0012285
/// HP:0000118["HP:0000118<br>
/// Phenotypic abnormality"]
/// HP:0000118 --> HP:0000707
/// HP:0000118 --> HP:0000818
/// HP:0000118 --> HP:0001939
/// HP:0000118 -----> HP:0002715
/// HP:0002011["HP:0002011<br>
/// Morphological central nervous system abnormality"]
/// HP:0002011 --> HP:0012443
/// HP:0031797["HP:0031797<br>
/// Clinical course"]
/// HP:0031797 --> HP:0003674
/// HP:0012639["HP:0012639<br>
/// Abnormal nervous system morphology"]
/// HP:0012639 --> HP:0002011
/// HP:0002715["HP:0002715<br>
/// Abnormality of the immune system"]
/// HP:0002715 --> HP:0010978
/// HP:0025454["HP:0025454<br>
/// Abnormal CSF metabolite concentration"]
/// HP:0000007["HP:0000007<br>
/// Autosomal recessive inheritance"]
/// HP:0000864["HP:0000864<br>
/// Abnormality of the hypothalamus-pituitary axis"]
/// HP:0000864 ---> HP:0012285
/// ```
#[derive(Default)]
pub struct Ontology {
    hpo_terms: Arena,
    genes: HashMap<GeneId, Gene>,
    omim_diseases: HashMap<OmimDiseaseId, OmimDisease>,
    hpo_version: (u16, u8, u8),
    categories: HpoGroup,
    modifier: HpoGroup,
}

impl Debug for Ontology {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Ontology with {} terns", self.hpo_terms.len())
    }
}

/// Public API of the Ontology
///
/// Those methods are all safe to use
impl Ontology {
    /// Initialize the [`Ontology`] from data provided by [Jax HPO](https://hpo.jax.org/)
    ///
    /// You must download:
    ///
    /// - Actual OBO data: [`hp.obo`](https://hpo.jax.org/app/data/ontology)
    /// - Links between HPO and OMIM diseases: [`phenotype.hpoa`](https://hpo.jax.org/app/data/annotations)
    /// - Links between HPO and Genes: [`phenotype_to_genes.txt`](http://purl.obolibrary.org/obo/hp/hpoa/phenotype_to_genes.txt)
    ///
    /// and then specify the folder where the data is stored.
    ///
    /// # Errors
    ///
    /// This method can fail for various reasons:
    ///
    /// - obo file not present or available: [`HpoError::CannotOpenFile`]
    /// - [`Ontology::add_gene`] failed
    /// - [`Ontology::add_omim_disease`] failed
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use hpo::Ontology;
    /// use hpo::HpoTermId;
    ///
    /// let ontology = Ontology::from_binary("/path/to/jax_hpo_data/").unwrap();
    ///
    /// assert!(ontology.len() == 26);
    ///
    /// let absent_term = HpoTermId::try_from("HP:9999999").unwrap();
    /// assert!(ontology.hpo(absent_term).is_none());
    ///
    /// let present_term = HpoTermId::try_from("HP:0000001").unwrap();
    /// let root_term = ontology.hpo(present_term).unwrap();
    /// assert_eq!(root_term.name(), "All");
    /// ```
    ///
    pub fn from_standard(folder: &str) -> HpoResult<Self> {
        let mut ont = Ontology::default();
        let path = Path::new(folder);
        let obo = path.join(crate::OBO_FILENAME);
        let gene = path.join(crate::GENE_FILENAME);
        let disease = path.join(crate::DISEASE_FILENAME);
        parser::load_from_standard_files(&obo, &gene, &disease, &mut ont)?;
        ont.calculate_information_content()?;
        ont.set_default_categories()?;
        ont.set_default_modifier()?;
        Ok(ont)
    }

    /// Build an Ontology from a binary data blob
    ///
    /// The data must be in the proper format, as defined in
    /// [`Ontology::as_bytes`]. This method adds all terms, creates the
    /// parent-child structure of the ontology, adds genes and Omim diseases
    /// and ensures proper inheritance of gene/disease annotations.
    /// It also calculates the `InformationContent` for every term.
    ///
    /// # Errors
    ///
    /// This method can fail for various reasons:
    ///
    /// - Binary file not available: [`HpoError::CannotOpenFile`]
    /// - `Ontology::add_genes_from_bytes` failed (TODO)
    /// - `Ontology::add_omim_disease_from_bytes` failed (TODO)
    /// - `add_terms_from_bytes` failed (TODO)
    /// - `add_parent_from_bytes` failed (TODO)
    /// - Size of binary data does not match the content: [`HpoError::ParseBinaryError`]
    ///
    /// # Examples
    ///
    /// ```
    /// use hpo::{Ontology, HpoTermId};
    ///
    /// let ontology = Ontology::from_binary("tests/example.hpo").unwrap();
    ///
    /// assert_eq!(ontology.len(), 26);
    ///
    /// let absent_term = HpoTermId::try_from("HP:9999999").unwrap();
    /// assert!(ontology.hpo(absent_term).is_none());
    ///
    /// let present_term = HpoTermId::try_from("HP:0000001").unwrap();
    /// let root_term = ontology.hpo(present_term).unwrap();
    /// assert_eq!(root_term.name(), "All");
    /// ```
    pub fn from_binary<P: AsRef<Path>>(filename: P) -> HpoResult<Self> {
        let bytes = match File::open(filename) {
            Ok(mut file) => {
                let len = file
                    .metadata()
                    .map_err(|_| {
                        HpoError::CannotOpenFile(
                            "unable to get filesize of binary file".to_string(),
                        )
                    })?
                    .len();
                let mut bytes = Vec::with_capacity(len.try_into()?);
                file.read_to_end(&mut bytes).map_err(|_| {
                    HpoError::CannotOpenFile("unable to read from binary file".to_string())
                })?;
                bytes
            }
            Err(_) => {
                return Err(crate::HpoError::CannotOpenFile(
                    "unable to open binary file".to_string(),
                ))
            }
        };
        Self::from_bytes(&bytes)
    }

    /// Build an Ontology from bytes
    ///
    /// The data must be in the proper format, as defined in
    /// [`Ontology::as_bytes`]. This method adds all terms, creates the
    /// parent-child structure of the ontology, adds genes and Omim diseases
    /// and ensures proper inheritance of gene/disease annotations.
    /// It also calculates the `InformationContent` for every term.
    ///
    /// # Errors
    ///
    /// This method can fail for various reasons:
    ///
    /// - Too few bytes or an invalid version
    /// - `Ontology::hpo_version_from_bytes` failed
    /// - `Ontology::add_genes_from_bytes` failed
    /// - `Ontology::add_omim_disease_from_bytes` failed
    /// - `add_terms_from_bytes` failed
    /// - `add_parent_from_bytes` failed
    /// - Size of binary data does not match the content: [`HpoError::ParseBinaryError`]
    ///
    ///
    /// # Examples
    ///
    /// ```
    /// use std::fs::File;
    /// use std::io::Read;
    /// use hpo::{Ontology, HpoTermId};
    ///
    /// let mut bytes = Vec::new();
    /// let mut file = File::open("tests/example.hpo").unwrap();
    /// file.read_to_end(&mut bytes).unwrap();
    /// let ontology = Ontology::from_bytes(&bytes).unwrap();
    ///
    /// assert_eq!(ontology.len(), 26);
    ///
    /// let absent_term = HpoTermId::try_from("HP:9999999").unwrap();
    /// assert!(ontology.hpo(absent_term).is_none());
    ///
    /// let present_term = HpoTermId::try_from("HP:0000001").unwrap();
    /// let root_term = ontology.hpo(present_term).unwrap();
    /// assert_eq!(root_term.name(), "All");
    /// ```
    pub fn from_bytes(bytes: &[u8]) -> HpoResult<Self> {
        let bytes = parser::binary::ontology::version(bytes)?;
        debug!("Parsing from bytes v{}", bytes.version());
        let mut ont = Ontology::default();

        let offset = ont.hpo_version_from_bytes(&bytes)?;

        let mut section_start = offset;
        let mut section_end: usize;

        // Terms
        let mut section_len = u32_from_bytes(&bytes[section_start..]) as usize;
        section_end = section_start + 4 + section_len;
        ont.add_terms_from_bytes(bytes.subset(section_start + 4..section_end));
        section_start += section_len + 4;

        // Term - Parents
        section_len = u32_from_bytes(&bytes[section_start..]) as usize;
        section_end += 4 + section_len;
        ont.add_parent_from_bytes(&bytes[section_start + 4..section_end]);
        ont.create_cache();
        section_start += section_len + 4;

        // Genes
        section_len = u32_from_bytes(&bytes[section_start..]) as usize;
        section_end += 4 + section_len;
        ont.add_genes_from_bytes(&bytes[section_start + 4..section_end])?;
        section_start += section_len + 4;

        // Omim Diseases
        section_len = u32_from_bytes(&bytes[section_start..]) as usize;
        section_end += 4 + section_len;
        ont.add_omim_disease_from_bytes(&bytes[section_start + 4..section_end])?;
        section_start += section_len + 4;

        if section_start == bytes.len() {
            ont.calculate_information_content()?;
            ont.set_default_categories()?;
            ont.set_default_modifier()?;
            Ok(ont)
        } else {
            Err(HpoError::ParseBinaryError)
        }
    }

    /// Returns a binary representation of the Ontology
    ///
    /// The binary data is separated into sections:
    ///
    /// - Metadata (HPO and Bindat Version) (see `Ontology::metadata_as_bytes`)
    /// - Terms (Names + IDs) (see `HpoTermInternal::as_bytes`)
    /// - Term - Parent connection (Child ID - Parent ID)
    ///   (see `HpoTermInternal::parents_as_byte`)
    /// - Genes (Names + IDs + Connected HPO Terms) ([`Gene::as_bytes`])
    /// - OMIM Diseases (Names + IDs + Connected HPO Terms)
    ///   ([`OmimDisease::as_bytes`])
    ///
    /// Every section starts with 4 bytes to indicate its size
    /// (big-endian encoded `u32`)
    ///
    /// This method is only useful if you use are modifying the ontology
    /// and want to save data for later re-use.
    ///
    /// # Panics
    ///
    /// Panics when the buffer length of any subsegment larger than `u32::MAX`
    ///
    /// # Examples
    ///
    /// ```
    /// use hpo::Ontology;
    /// let ontology = Ontology::from_binary("tests/example.hpo").unwrap();
    /// let bytes = ontology.as_bytes();
    /// ```
    pub fn as_bytes(&self) -> Vec<u8> {
        fn usize_to_u32(n: usize) -> u32 {
            n.try_into().expect("unable to convert {n} to u32")
        }
        let mut res = Vec::new();

        // Add metadata, version info
        res.append(&mut self.metadata_as_bytes());

        // All HPO Terms
        let mut buffer = Vec::new();
        for term in self.hpo_terms.values() {
            buffer.append(&mut term.as_bytes());
        }
        res.append(&mut usize_to_u32(buffer.len()).to_be_bytes().to_vec());
        res.append(&mut buffer);

        // All Term - Parent connections
        buffer.clear();
        for term in self.hpo_terms.values() {
            buffer.append(&mut term.parents_as_byte());
        }
        res.append(&mut usize_to_u32(buffer.len()).to_be_bytes().to_vec());
        res.append(&mut buffer);

        // Genes and Gene-Term connections
        buffer.clear();
        for gene in self.genes.values() {
            buffer.append(&mut gene.as_bytes());
        }
        res.append(&mut usize_to_u32(buffer.len()).to_be_bytes().to_vec());
        res.append(&mut buffer);

        // OMIM Disease and Disease-Term connections
        buffer.clear();
        for omim_disease in self.omim_diseases.values() {
            buffer.append(&mut omim_disease.as_bytes());
        }
        res.append(&mut usize_to_u32(buffer.len()).to_be_bytes().to_vec());
        res.append(&mut buffer);

        res
    }

    /// Returns the number of HPO-Terms in the Ontology
    ///
    /// # Examples
    ///
    /// ```
    /// use hpo::Ontology;
    /// let ontology = Ontology::from_binary("tests/example.hpo").unwrap();
    /// assert_eq!(ontology.len(), 26);
    /// ```
    pub fn len(&self) -> usize {
        self.hpo_terms.len()
    }

    /// Returns `true` if the Ontology does not contain any HPO-Terms
    ///
    /// # Examples
    ///
    /// ```
    /// use hpo::Ontology;
    /// let ontology = Ontology::default();
    /// assert!(ontology.is_empty());
    /// ```
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Returns the [`HpoTerm`] of the provided [`HpoTermId`]
    ///
    /// If no such term is present in the Ontolgy, `None` is returned
    ///
    /// # Examples
    ///
    /// ```
    /// use hpo::Ontology;
    /// let ontology = Ontology::from_binary("tests/example.hpo").unwrap();
    /// let term = ontology.hpo(11017u32).unwrap();
    /// assert_eq!(term.name(), "Abnormal cellular physiology");
    /// assert!(ontology.hpo(66666u32).is_none());
    /// ```
    pub fn hpo<I: Into<HpoTermId>>(&self, term_id: I) -> Option<HpoTerm> {
        HpoTerm::try_new(self, term_id.into()).ok()
    }

    /// Returns an Iterator of all [`HpoTerm`]s from the Ontology
    ///
    /// # Examples
    ///
    /// ```
    /// use hpo::Ontology;
    /// let ontology = Ontology::from_binary("tests/example.hpo").unwrap();
    /// for term in ontology.hpos() {
    ///     println!("{}", term.name());
    /// }
    /// ```
    ///
    pub fn hpos(&self) -> Iter<'_> {
        self.into_iter()
    }

    /// Returns a reference to the [`Gene`] of the provided [`GeneId`]
    ///
    /// If no such gene is present, `None` is returned
    ///
    /// # Examples
    ///
    /// ```
    /// use hpo::Ontology;
    /// let ontology = Ontology::from_binary("tests/example.hpo").unwrap();
    /// let gene = ontology.gene(&57505u32.into()).unwrap();
    /// assert_eq!(gene.name(), "AARS2");
    /// ```
    pub fn gene(&self, gene_id: &GeneId) -> Option<&Gene> {
        self.genes.get(gene_id)
    }

    /// Returns a reference to the [`Gene`] with the provided symbol / name
    ///
    /// If no such gene is present, `None` is returned
    ///
    /// # Note
    ///
    /// `Gene`s are not index by name, so this method searches through all
    /// genes. If you can, prefer using [`Ontology::gene`] with the [`GeneId`].
    ///
    /// # Examples
    ///
    /// ```
    /// use hpo::Ontology;
    /// let ontology = Ontology::from_binary("tests/example.hpo").unwrap();
    ///
    /// let gene = ontology.gene_by_name("AARS2").unwrap();
    /// assert_eq!(gene.name(), "AARS2");
    ///
    /// assert!(ontology.gene_by_name("FOOBAR66").is_none());
    /// ```
    pub fn gene_by_name(&self, symbol: &str) -> Option<&Gene> {
        self.genes.values().find(|&gene| gene.name() == symbol)
    }

    /// Returns an Iterator of all [`Gene`]s from the Ontology
    ///
    /// It is likely that the return type will change to a dedicated Iterator
    ///
    /// # Examples
    ///
    /// ```
    /// use hpo::Ontology;
    /// let ontology = Ontology::from_binary("tests/example.hpo").unwrap();
    /// for gene in ontology.genes() {
    ///     println!("{}", gene.name());
    /// }
    /// ```
    pub fn genes(&self) -> std::collections::hash_map::Values<'_, GeneId, Gene> {
        self.genes.values()
    }

    /// Returns a reference to the [`OmimDisease`] of the provided [`OmimDiseaseId`]
    ///
    /// If no such disease is present, `None` is returned
    ///
    /// # Examples
    ///
    /// ```
    /// use hpo::Ontology;
    /// let ontology = Ontology::from_binary("tests/example.hpo").unwrap();
    /// let disease = ontology.omim_disease(&601495u32.into()).unwrap();
    /// assert_eq!(disease.name(), "Agammaglobulinemia 1, autosomal recessive");
    /// ```
    pub fn omim_disease(&self, omim_disease_id: &OmimDiseaseId) -> Option<&OmimDisease> {
        self.omim_diseases.get(omim_disease_id)
    }

    /// Returns an Iterator of all [`OmimDisease`]s from the Ontology
    ///
    /// It is likely that the return type will change to a dedicated Iterator
    ///
    /// # Examples
    ///
    /// ```
    /// use hpo::Ontology;
    /// let ontology = Ontology::from_binary("tests/example.hpo").unwrap();
    /// for disease in ontology.omim_diseases() {
    ///     println!("{}", disease.name());
    /// }
    /// ```
    pub fn omim_diseases(
        &self,
    ) -> std::collections::hash_map::Values<'_, OmimDiseaseId, OmimDisease> {
        self.omim_diseases.values()
    }

    /// Returns the Jax-Ontology release version
    ///
    /// e.g. `2023-03-13`
    pub fn hpo_version(&self) -> String {
        format!(
            "{:0>4}-{:0>2}-{:0>2}",
            self.hpo_version.0, self.hpo_version.1, self.hpo_version.2,
        )
    }

    /// Compares `self` to another `Ontology` to identify added/removed terms, genes and diseases
    ///
    /// # Examples
    ///
    /// ```
    /// use hpo::Ontology;
    ///
    /// let ontology_1 = Ontology::from_binary("tests/example.hpo").unwrap();
    /// let mut ontology_2 = Ontology::default();
    ///
    /// ontology_2.add_gene("FOOBAR", "666666").unwrap();
    ///
    /// let compare = ontology_1.compare(&ontology_2);
    /// assert_eq!(compare.added_hpo_terms().len(), 0);
    /// assert_eq!(compare.removed_hpo_terms().len(), 26);
    /// assert_eq!(compare.added_genes().len(), 1);
    /// ```
    pub fn compare<'a>(&'a self, other: &'a Ontology) -> Comparison {
        Comparison::new(self, other)
    }

    /// Constructs a smaller ontology that contains only the `leaves` terms and
    /// all terms needed to connect to each leaf to `root`
    ///
    /// # Errors
    ///
    /// Fails if `root` is not an ancestor of all leaves
    ///
    /// # Examples
    ///
    /// ```
    /// use hpo::Ontology;
    ///
    /// let ontology = Ontology::from_binary("tests/example.hpo").unwrap();
    /// let ontology_2 = ontology.sub_ontology(
    ///     ontology.hpo(118u32).unwrap(),
    ///     vec![ontology.hpo(11017u32).unwrap()]
    /// ).unwrap();
    ///
    /// assert_eq!(ontology_2.len(), 3);
    /// ```
    pub fn sub_ontology<'a, T: IntoIterator<Item = HpoTerm<'a>>>(
        &self,
        root: HpoTerm,
        leaves: T,
    ) -> Result<Self, HpoError> {
        let mut terms = HashSet::new();
        for term in leaves {
            terms.insert(self.get_unchecked(term.id()));
            for parent in term
                .path_to_ancestor(&root)
                .ok_or(HpoError::NotImplemented)?
            {
                terms.insert(self.get_unchecked(parent));
            }
        }
        let ids: HpoGroup = terms.iter().map(|term| *term.id()).collect();

        let mut ont = Self::default();
        for &term in &terms {
            let mut copied_term = HpoTermInternal::new(term.name().to_string(), *term.id());
            *copied_term.obsolete_mut() = term.obsolete();
            *copied_term.replacement_mut() = term.replacement();
            ont.add_term(copied_term);
        }
        for term in &terms {
            for parent in term.parents() {
                if ids.contains(&parent) {
                    ont.add_parent(parent, *term.id());
                }
            }
        }

        ont.create_cache();

        // Iterate all genes, check the associated terms and see if one
        // is part of the `ids` set
        for gene in self.genes() {
            let matched_terms = gene.hpo_terms() & &ids;
            if matched_terms.is_empty() {
                continue;
            }
            let gene_id = ont.add_gene(
                self.gene(gene.id()).ok_or(HpoError::DoesNotExist)?.name(),
                &gene.id().as_u32().to_string(),
            )?;
            for term in &matched_terms {
                ont.link_gene_term(term, gene_id)?;
                ont.gene_mut(&gene_id)
                    .ok_or(HpoError::DoesNotExist)?
                    .add_term(term);
            }
        }

        // Iterate all genes, check the associated terms and see if one
        // is part of the `ids` set
        for omim_disease in self.omim_diseases() {
            let matched_terms = omim_disease.hpo_terms() & &ids;
            if matched_terms.is_empty() {
                continue;
            }
            let omim_disease_id = ont.add_omim_disease(
                self.omim_disease(omim_disease.id())
                    .ok_or(HpoError::DoesNotExist)?
                    .name(),
                &omim_disease.id().as_u32().to_string(),
            )?;
            for term in &matched_terms {
                ont.link_omim_disease_term(term, omim_disease_id)?;
                ont.omim_disease_mut(&omim_disease_id)
                    .ok_or(HpoError::DoesNotExist)?
                    .add_term(term);
            }
        }
        ont.calculate_information_content()?;

        Ok(ont)
    }

    /// Returns the code to crate a `Mermaid` flow diagram
    ///
    /// This is meant to be used with smaller ontologies, e.g. from [`Ontology::sub_ontology`]
    pub fn as_mermaid(&self) -> String {
        let mut code = String::new();
        code.push_str("graph TD\n");
        for term in self {
            code.push_str(&format!(
                "{}[\"{}\n{}\"]\n",
                term.id(),
                term.id(),
                term.name()
            ));
            for child in term.children() {
                code.push_str(&format!("{} --> {}\n", term.id(), child.id()));
            }
        }
        code
    }

    /// Returns a reference to the categories of the Ontology
    ///
    /// Categories are top-level `HpoTermId`s used for
    /// categorizing individual `HpoTerm`s.
    ///
    /// See [`Ontology::set_default_categories()`] for more information
    ///
    /// # Examples
    ///
    /// ```
    /// use hpo::{HpoTerm, Ontology};
    ///
    /// let mut ontology = Ontology::from_binary("tests/example.hpo").unwrap();
    /// assert_eq!(ontology.categories().len(), 6);
    /// ```
    pub fn categories(&self) -> &HpoGroup {
        &self.categories
    }

    /// Returns a mutable reference to the categories vector
    ///
    /// This is a vector that should contain top-level `HpoTermId`s used for
    /// categorizing individual `HpoTerm`s.
    ///
    /// See [`Ontology::set_default_categories()`]
    pub fn categories_mut(&mut self) -> &mut HpoGroup {
        &mut self.categories
    }

    /// Sets the default categories for the Ontology
    ///
    /// By default, each direct child of [`Phenotypic abnormality`](crate::PHENOTYPE_ID)
    /// is considered one category, e.g.:
    ///
    /// - `HP:0000152 | Abnormality of head or neck`
    /// - `HP:0001507 | Growth abnormality`
    /// - ...
    ///
    /// In additon to all other direct children of `HP:0000001 | All`, e.g.:
    ///
    /// - `HP:0000005 | Mode of inheritance`
    /// - `HP:0012823 | Clinical modifier`
    /// - ...
    ///
    /// # Errors
    ///
    /// This method requires that the main-category terms:
    ///
    /// - `HP:0000001 | All`
    /// - `HP:0000118 | Phenotypic abnormality`
    ///
    /// are present in the Ontology.
    pub fn set_default_categories(&mut self) -> HpoResult<()> {
        let root = self.hpo(1u32.into()).ok_or(HpoError::DoesNotExist)?;
        let phenotypes = self
            .hpo(crate::PHENOTYPE_ID)
            .ok_or(HpoError::DoesNotExist)?;
        self.categories = root
            .children_ids()
            .iter()
            .filter(|id| id != &crate::PHENOTYPE_ID)
            .chain(phenotypes.children_ids())
            .collect();
        Ok(())
    }

    /// Returns a reference to the modifier root terms of the Ontology
    ///
    /// See [`Ontology::set_default_modifier()`] for more information
    ///
    /// # Examples
    ///
    /// ```
    /// use hpo::{HpoTerm, Ontology};
    ///
    /// let mut ontology = Ontology::from_binary("tests/example.hpo").unwrap();
    /// assert_eq!(ontology.modifier().len(), 2);
    /// ```
    pub fn modifier(&self) -> &HpoGroup {
        &self.modifier
    }

    /// Returns a mutable reference to the modifier vector
    ///
    /// This is a vector that should contain top-level `HpoTermId`s that are
    /// representing modifier terms.
    ///
    /// See [`Ontology::set_default_modifier()`]
    pub fn modifier_mut(&mut self) -> &mut HpoGroup {
        &mut self.modifier
    }

    /// Sets the default modifier categories for the Ontology
    ///
    /// The default is very opinionated and declares everything a modifier
    /// that is not part of the [`Phenotypic abnormality`](crate::PHENOTYPE_ID)
    /// category.
    ///
    /// # Errors
    ///
    /// This method requires that the root term `HP:0000001 | All` is present.
    pub fn set_default_modifier(&mut self) -> HpoResult<()> {
        self.modifier = self
            .hpo(1u32.into())
            .ok_or(HpoError::DoesNotExist)?
            .children_ids()
            .iter()
            .filter(|id| id != &crate::PHENOTYPE_ID)
            .collect();
        Ok(())
    }
}

/// Methods to add annotations
///
/// These methods should rarely (if ever) be used by clients.
/// Calling these functions might disrupt the Ontology and associated terms.
impl Ontology {
    /// Crates and inserts a new term to the ontology
    ///
    /// This method does not link the term to its parents or to any annotations
    ///
    /// # Examples
    ///
    /// ```
    /// use hpo::Ontology;
    ///
    /// let mut ontology = Ontology::default();
    /// ontology.insert_term("FooBar".into(), 1u32);
    ///
    /// assert_eq!(ontology.len(), 1);
    /// ```
    pub fn insert_term<I: Into<HpoTermId>>(&mut self, name: String, id: I) {
        let term = HpoTermInternal::new(name, id.into());
        self.hpo_terms.insert(term);
    }

    /// Add a connection from an [`HpoTerm`] to its parent
    ///
    /// This method is called once for every dependency in the Ontology during the initialization.
    ///
    /// There should rarely be a need to call this method outside of the ontology building
    ///
    /// # Panics
    ///
    /// This method will panic if the `parent_id` or `child_id` is not present in the Ontology
    ///
    /// # Examples
    ///
    /// ```
    /// use hpo::Ontology;
    ///
    /// let mut ontology = Ontology::default();
    /// ontology.insert_term("Foo".into(), 1u32);
    /// ontology.insert_term("Bar".into(), 2u32);
    ///
    /// ontology.add_parent(1u32, 2u32);
    ///
    /// assert!(ontology.hpo(2u32).unwrap().parent_ids().contains(&1u32.into()));
    /// ```
    pub fn add_parent<I: Into<HpoTermId> + Copy, J: Into<HpoTermId> + Copy>(
        &mut self,
        parent_id: I,
        child_id: J,
    ) {
        let parent = self.get_unchecked_mut(parent_id);
        parent.add_child(child_id);

        let child = self.get_unchecked_mut(child_id);
        child.add_parent(parent_id);
    }

    /// Crates and caches the `all_parents` values for every term
    ///
    /// This method can only be called once and afterwards no new terms
    /// should be added to the Ontology anymore and no new term-parent connection
    /// should be created.
    /// Since this method caches the results, rerunning it will not cause a new
    /// calculation.
    ///
    /// # Examples
    ///
    /// ```
    /// use hpo::Ontology;
    ///
    /// let mut ontology = Ontology::default();
    /// ontology.insert_term("Root".into(), 1u32);
    /// ontology.insert_term("Foo".into(), 2u32);
    /// ontology.insert_term("Bar".into(), 3u32);
    ///
    /// ontology.add_parent(1u32, 2u32);
    /// ontology.add_parent(2u32, 3u32);
    ///
    /// // At this point #3 does not have info about grandparents
    /// assert!(!ontology.hpo(3u32).unwrap().all_parent_ids().contains(&1u32.into()));
    ///
    /// ontology.create_cache();
    /// assert!(ontology.hpo(3u32).unwrap().all_parent_ids().contains(&1u32.into()));
    /// ```
    pub fn create_cache(&mut self) {
        let term_ids: Vec<HpoTermId> = self.hpo_terms.keys();

        for id in term_ids {
            self.create_cache_of_grandparents(id);
        }
    }

    /// Add a gene to the Ontology. and return the [`GeneId`]
    ///
    /// If the gene does not yet exist, a new [`Gene`] entity is created
    /// and stored in the Ontology.
    /// If the gene already exists in the ontology, it is not added again.
    ///
    /// # Note
    ///
    /// Adding a gene does not connect it to any HPO terms.
    /// Use [`Ontology::link_gene_term`] for creating connections.
    ///
    /// # Errors
    ///
    /// If the `gene_id` is invalid, an [`HpoError::ParseIntError`] is returned
    ///
    /// # Examples
    ///
    /// ```
    /// use hpo::Ontology;
    ///
    /// let mut ontology = Ontology::default();
    /// assert!(ontology.gene(&1u32.into()).is_none());
    ///
    /// ontology.add_gene("Foo", "1");
    ///
    /// // Genes can be iterated...
    /// let mut gene_iterator = ontology.genes();
    /// let gene = gene_iterator.next().unwrap();
    /// assert_eq!(gene.name(), "Foo");
    /// assert!(gene_iterator.next().is_none());
    ///
    /// // .. or accessed directly
    /// assert!(ontology.gene(&1u32.into()).is_some());
    /// ```
    pub fn add_gene(&mut self, gene_name: &str, gene_id: &str) -> HpoResult<GeneId> {
        let id = GeneId::try_from(gene_id)?;
        match self.genes.entry(id) {
            std::collections::hash_map::Entry::Occupied(_) => Ok(id),
            std::collections::hash_map::Entry::Vacant(entry) => {
                entry.insert(Gene::new(id, gene_name));
                Ok(id)
            }
        }
    }

    /// Add a OMIM disease to the Ontology. and return the [`OmimDiseaseId`]
    ///
    /// If the disease does not yet exist, a new [`OmimDisease`] entity is
    /// created and stored in the Ontology.
    /// If the disease already exists in the ontology, it is not added again.
    ///
    /// # Note
    ///
    /// Adding a disease does not connect it to any HPO terms.
    /// Use [`Ontology::link_omim_disease_term`] for creating connections.
    ///
    /// # Errors
    ///
    /// If the `omim_disease_id` is invalid, an [`HpoError::ParseIntError`] is returned
    ///
    /// # Examples
    ///
    /// ```
    /// use hpo::Ontology;
    ///
    /// let mut ontology = Ontology::default();
    /// assert!(ontology.omim_disease(&1u32.into()).is_none());
    ///
    /// ontology.add_omim_disease("Foo", "1");
    ///
    /// // Diseases can be iterated...
    /// let mut disease_iterator = ontology.omim_diseases();
    /// let omim_disease = disease_iterator.next().unwrap();
    /// assert_eq!(omim_disease.name(), "Foo");
    /// assert!(disease_iterator.next().is_none());
    ///
    /// // .. or accessed directly
    /// assert!(ontology.omim_disease(&1u32.into()).is_some());
    /// ```
    pub fn add_omim_disease(
        &mut self,
        omim_disease_name: &str,
        omim_disease_id: &str,
    ) -> HpoResult<OmimDiseaseId> {
        let id = OmimDiseaseId::try_from(omim_disease_id)?;
        match self.omim_diseases.entry(id) {
            std::collections::hash_map::Entry::Occupied(_) => Ok(id),
            std::collections::hash_map::Entry::Vacant(entry) => {
                entry.insert(OmimDisease::new(id, omim_disease_name));
                Ok(id)
            }
        }
    }

    /// Add the [`Gene`] as annotation to the [`HpoTerm`]
    ///
    /// The gene will be recursively connected to all parent `HpoTerms` as well.
    ///
    /// This method does not add the HPO-term to the [`Gene`], this must be handled
    /// by the client.
    ///
    /// # Errors
    ///
    /// If the HPO term is not present, an [`HpoError::DoesNotExist`] is returned
    ///
    /// # Examples
    ///
    /// ```
    /// use hpo::Ontology;
    /// use hpo::annotations::GeneId;
    ///
    /// let mut ontology = Ontology::default();
    /// ontology.insert_term("Term-Foo".into(), 1u32);
    /// ontology.add_gene("Foo", "5");
    /// ontology.link_gene_term(1u32, GeneId::from(5u32)).unwrap();
    ///
    /// let term = ontology.hpo(1u32).unwrap();
    /// assert_eq!(term.genes().next().unwrap().name(), "Foo");
    /// ```
    pub fn link_gene_term<I: Into<HpoTermId>>(
        &mut self,
        term_id: I,
        gene_id: GeneId,
    ) -> HpoResult<()> {
        let term = self.get_mut(term_id).ok_or(HpoError::DoesNotExist)?;

        if term.add_gene(gene_id) {
            // If the gene is already associated to the term, this branch will
            // be skipped. That is desired, because by definition
            // all parent terms are already linked as well
            let parents = term.all_parents().clone();
            for parent in &parents {
                self.link_gene_term(parent, gene_id)?;
            }
        }
        Ok(())
    }

    /// Add the [`OmimDisease`] as annotation to the [`HpoTerm`]
    ///
    /// The disease will be recursively connected to all parent `HpoTerms` as well.
    ///
    /// This method does not add the HPO-term to the [`OmimDisease`], this
    /// must be handled by the client.
    ///
    /// # Errors
    ///
    /// If the HPO term is not present, an [`HpoError`] is returned
    ///
    /// # Examples
    ///
    /// ```
    /// use hpo::Ontology;
    /// use hpo::annotations::OmimDiseaseId;
    ///
    /// let mut ontology = Ontology::default();
    /// ontology.insert_term("Term-Foo".into(), 1u32);
    /// ontology.add_omim_disease("Foo", "5");
    /// ontology.link_omim_disease_term(1u32, OmimDiseaseId::from(5u32)).unwrap();
    ///
    /// let term = ontology.hpo(1u32).unwrap();
    /// assert_eq!(term.omim_diseases().next().unwrap().name(), "Foo");
    /// ```
    pub fn link_omim_disease_term<I: Into<HpoTermId>>(
        &mut self,
        term_id: I,
        omim_disease_id: OmimDiseaseId,
    ) -> HpoResult<()> {
        let term = self.get_mut(term_id).ok_or(HpoError::DoesNotExist)?;

        if term.add_omim_disease(omim_disease_id) {
            // If the disease is already associated to the term, this branch will
            // be skipped. That is desired, because by definition
            // all parent terms are already linked as well
            let parents = term.all_parents().clone();
            for parent in &parents {
                self.link_omim_disease_term(parent, omim_disease_id)?;
            }
        }
        Ok(())
    }

    /// Returns a mutable reference to the [`Gene`] of the provided [`GeneId`]
    ///
    /// If no such gene is present, `None` is returned
    ///
    /// # Examples
    ///
    /// ```
    /// use hpo::Ontology;
    ///
    /// let mut ontology = Ontology::from_binary("tests/example.hpo").unwrap();
    ///
    /// let mut gene = ontology.gene_mut(&57505u32.into()).unwrap();
    /// assert_eq!(gene.hpo_terms().len(), 10);
    /// gene.add_term(1u32);
    /// assert_eq!(gene.hpo_terms().len(), 11);
    /// ```
    pub fn gene_mut(&mut self, gene_id: &GeneId) -> Option<&mut Gene> {
        self.genes.get_mut(gene_id)
    }

    /// Returns a mutable reference to the [`OmimDisease`] of the provided [`OmimDiseaseId`]
    ///
    /// If no such disease is present, `None` is returned
    ///
    /// # Examples
    ///
    /// ```
    /// use hpo::Ontology;
    ///
    /// let mut ontology = Ontology::from_binary("tests/example.hpo").unwrap();
    ///
    /// let mut disease = ontology.omim_disease_mut(&601495u32.into()).unwrap();
    /// assert_eq!(disease.hpo_terms().len(), 1);
    /// disease.add_term(1u32);
    /// assert_eq!(disease.hpo_terms().len(), 2);
    /// ```
    pub fn omim_disease_mut(
        &mut self,
        omim_disease_id: &OmimDiseaseId,
    ) -> Option<&mut OmimDisease> {
        self.omim_diseases.get_mut(omim_disease_id)
    }

    /// Calculates the [`crate::term::InformationContent`]s for every term
    ///
    /// This method should only be called **after** all terms are added,
    /// connected and all genes and diseases are linked as well.
    ///
    /// It can be called repeatedly, all values are recalculated each time,
    /// as long as the Ontology contains at least 1 gene/disease.
    /// When no genes/diseases are present, the IC is not calculated nor updated.
    ///
    /// # Errors
    ///
    /// This method returns an error if there are more Genes or Terms than `u16::MAX`
    /// because larger numbers can't be safely converted to `f32`
    ///
    /// # Examples
    ///
    /// ```
    /// use hpo::Ontology;
    ///
    /// let mut ontology = Ontology::default();
    ///
    /// // [all kind of logic to add terms, diseases, genes....]
    ///
    /// ontology.calculate_information_content().unwrap();
    /// ```
    pub fn calculate_information_content(&mut self) -> HpoResult<()> {
        self.calculate_gene_ic()?;
        self.calculate_omim_disease_ic()?;
        Ok(())
    }
}

/// Crate-only functions for setting up and building the Ontology
///
/// Those methods should not be exposed publicly
impl Ontology {
    /// Insert an `HpoTermInternal` to the ontology
    ///
    /// This method does not link the term to its parents or to any annotations
    pub(crate) fn add_term(&mut self, term: HpoTermInternal) -> HpoTermId {
        let id = *term.id();
        self.hpo_terms.insert(term);
        id
    }

    pub(crate) fn set_hpo_version(&mut self, version: (u16, u8, u8)) {
        self.hpo_version = version;
    }

    /// Parses `Bytes` into the Jax-Ontology release version
    fn hpo_version_from_bytes(&mut self, bytes: &Bytes) -> HpoResult<usize> {
        if bytes.version() == BinaryVersion::V1 {
            self.set_hpo_version((0u16, 0u8, 0u8));
            Ok(0)
        } else {
            if bytes.len() < 4 {
                return Err(HpoError::ParseBinaryError);
            }
            let year = u16::from_be_bytes([bytes[0], bytes[1]]);
            let month = u8::from_be_bytes([bytes[2]]);
            let day = u8::from_be_bytes([bytes[3]]);
            self.set_hpo_version((year, month, day));
            Ok(4)
        }
    }

    /// Returns a binary representation of the Ontology's metadata
    ///
    /// It adds the HPO-identifying bytes `HPO`, the version
    /// and the Jax-Ontology release version
    fn metadata_as_bytes(&self) -> Vec<u8> {
        let mut bytes = Vec::new();
        // HPO header
        bytes.extend_from_slice(&[0x48, 0x50, 0x4f]);

        // Version
        bytes.push(0x2);

        bytes.extend_from_slice(&self.hpo_version.0.to_be_bytes()[..]);
        bytes.push(self.hpo_version.1);
        bytes.push(self.hpo_version.2);
        bytes
    }

    /// Adds an [`HpoTerm`] to the ontology
    ///
    /// This method is part of the Ontology-building, based on the binary
    /// data format and requires a specified data layout.
    ///
    /// The method assumes that the data is in the right format and also
    /// assumes that the caller takes care of handling all consistencies
    /// like parent-child connection etc.
    ///
    /// See [`HpoTermInternal::as_bytes`] for explanation of the binary layout.
    fn add_terms_from_bytes(&mut self, bytes: Bytes) {
        for term in BinaryTermBuilder::new(bytes) {
            self.add_term(term);
        }
    }

    /// Connects an [`HpoTerm`] to its parent term
    ///
    /// This method is part of the Ontology-building, based on the binary
    /// data format and requires a specified data layout.
    ///
    /// The method assumes that the data is in the right format and also
    /// assumes that the caller will populate the `all_parents` caches for
    /// each term.
    ///
    /// See [`HpoTermInternal::parents_as_byte`] for explanation of the binary layout.
    ///
    /// # Panics
    ///
    /// This method will panic if the length of bytes does not exactly correspond
    /// to the contained data
    fn add_parent_from_bytes(&mut self, bytes: &[u8]) {
        let mut idx: usize = 0;
        loop {
            if idx == bytes.len() {
                break;
            }
            let n_parents = u32_from_bytes(&bytes[idx..]) as usize;

            idx += 4;
            let term =
                HpoTermId::from([bytes[idx], bytes[idx + 1], bytes[idx + 2], bytes[idx + 3]]);
            idx += 4;
            for _ in 0..n_parents {
                let parent =
                    HpoTermId::from([bytes[idx], bytes[idx + 1], bytes[idx + 2], bytes[idx + 3]]);
                self.add_parent(parent, term);
                idx += 4;
            }
        }
    }

    /// Adds genes to the ontoloigy and connects them to connected terms
    ///
    /// This method is part of the Ontology-building, based on the binary
    /// data format and requires a specified data layout.
    ///
    /// It connects all connected terms and their parents properly. The
    /// method assumes that the bytes encode all gene-term connections.
    ///
    /// See [`Gene::as_bytes`] for explanation of the binary layout
    fn add_genes_from_bytes(&mut self, bytes: &[u8]) -> HpoResult<()> {
        let mut idx: usize = 0;
        loop {
            if idx >= bytes.len() {
                break;
            }
            let gene_len = u32_from_bytes(&bytes[idx..]) as usize;
            let gene = Gene::try_from(&bytes[idx..idx + gene_len])?;
            for term in gene.hpo_terms() {
                self.link_gene_term(term, *gene.id())?;
            }
            self.genes.insert(*gene.id(), gene);
            idx += gene_len;
        }
        Ok(())
    }

    /// Adds [`OmimDisease`]s to the ontoloigy and connects them to connected terms
    ///
    /// This method is part of the Ontology-building, based on the binary
    /// data format and requires a specified data layout.
    ///
    /// It connects all connected terms and their parents properly. The
    /// method assumes that the bytes encode all Disease-term connections.
    ///
    /// See [`OmimDisease::as_bytes`] for explanation of the binary layout
    fn add_omim_disease_from_bytes(&mut self, bytes: &[u8]) -> HpoResult<()> {
        let mut idx: usize = 0;
        loop {
            if idx >= bytes.len() {
                break;
            }
            let disease_len = u32_from_bytes(&bytes[idx..]) as usize;
            let disease = OmimDisease::try_from(&bytes[idx..idx + disease_len])?;
            for term in disease.hpo_terms() {
                self.link_omim_disease_term(term, *disease.id())?;
            }
            self.omim_diseases.insert(*disease.id(), disease);
            idx += disease_len;
        }
        Ok(())
    }

    /// This method is part of the cache creation to link all terms to their
    /// direct and indirect parents (grandparents)
    ///
    /// # Panics
    ///
    /// This method will panic if the `term_id` is not present in the Ontology
    fn all_grandparents(&mut self, term_id: HpoTermId) -> &HpoGroup {
        if !self.get_unchecked(term_id).parents_cached() {
            self.create_cache_of_grandparents(term_id);
        }
        let term = self.get_unchecked(term_id);
        term.all_parents()
    }

    /// This method is part of the cache creation to link all terms to their
    /// direct and indirect parents (grandparents)
    ///
    /// It will (somewhat) recursively iterate all parents and copy all their parents.
    /// During this recursion, the list of `all_parents` is cached in each term that was
    /// iterated.
    ///
    /// The logic is that the recursion bubbles up all the way to the top of the ontolgy
    /// and then caches the list of direct and indirect parents for every term bubbling
    /// back down. The recursion does not reach the top level again, because it will stop
    /// once it reaches a term with already cached `all_parents`.
    ///
    /// # Panics
    ///
    /// This method will panic if the `term_id` is not present in the Ontology
    fn create_cache_of_grandparents(&mut self, term_id: HpoTermId) {
        let mut res = HpoGroup::default();
        let parents = self.get_unchecked(term_id).parents().clone();
        for parent in &parents {
            let grandparents = self.all_grandparents(parent);
            for gp in grandparents {
                res.insert(gp);
            }
        }
        let term = self.get_unchecked_mut(term_id);
        *term.all_parents_mut() = res.bitor(&parents);
    }

    /// Returns the `HpoTermInternal` with the given `HpoTermId`
    ///
    /// Returns `None` if no such term is present
    pub(crate) fn get<I: Into<HpoTermId>>(&self, term_id: I) -> Option<&HpoTermInternal> {
        self.hpo_terms.get(term_id.into())
    }

    /// Returns the `HpoTermInternal` with the given `HpoTermId`
    ///
    /// This method should only be called if the caller is sure that the term actually
    /// exists, e.g. during an iteration of all `HpoTermId`s.
    ///
    /// # Panics
    ///
    /// This method will panic if the `term_id` is not present in the Ontology
    pub(crate) fn get_unchecked<I: Into<HpoTermId>>(&self, term_id: I) -> &HpoTermInternal {
        self.hpo_terms.get_unchecked(term_id.into())
    }

    /// Returns a mutable reference to the `HpoTermInternal` with the given `HpoTermId`
    ///
    /// Returns `None` if no such term is present
    fn get_mut<I: Into<HpoTermId>>(&mut self, term_id: I) -> Option<&mut HpoTermInternal> {
        self.hpo_terms.get_mut(term_id.into())
    }

    /// Returns a mutable reference to the `HpoTermInternal` with the given `HpoTermId`
    ///
    /// This method should only be called if the caller is sure that the term actually
    /// exists, e.g. during an iteration of all `HpoTermId`s.
    ///
    /// # Panics
    ///
    /// This method will panic if the `term_id` is not present in the Ontology
    fn get_unchecked_mut<I: Into<HpoTermId>>(&mut self, term_id: I) -> &mut HpoTermInternal {
        self.hpo_terms.get_unchecked_mut(term_id.into())
    }

    /// Calculates the gene-specific Information Content for every term
    ///
    /// If no genes are present in the Ontology, no IC are calculated
    fn calculate_gene_ic(&mut self) -> HpoResult<()> {
        let n_genes = self.genes.len();
        for term in self.hpo_terms.values_mut() {
            let current_genes = term.genes().len();
            term.information_content_mut()
                .set_gene(n_genes, current_genes)?;
        }
        Ok(())
    }

    /// Calculates the Omim-Disease-specific Information Content for every term
    ///
    /// If no diseases are present in the Ontology, no IC are calculated
    fn calculate_omim_disease_ic(&mut self) -> HpoResult<()> {
        let n_omim_diseases = self.omim_diseases.len();

        for term in self.hpo_terms.values_mut() {
            let current_diseases = term.omim_diseases().len();
            term.information_content_mut()
                .set_omim_disease(n_omim_diseases, current_diseases)?;
        }
        Ok(())
    }
}

/// Iterates the Ontology and yields [`HpoTerm`]s
pub struct Iter<'a> {
    inner: termarena::Iter<'a>,
    ontology: &'a Ontology,
}

impl<'a> std::iter::Iterator for Iter<'a> {
    type Item = HpoTerm<'a>;
    fn next(&mut self) -> Option<Self::Item> {
        match self.inner.next() {
            Some(term) => Some(
                HpoTerm::try_new(self.ontology, term)
                    .expect("Iterator can only iterate valid HpoTermIds"),
            ),
            None => None,
        }
    }
}

impl<'a> IntoIterator for &'a Ontology {
    type Item = HpoTerm<'a>;
    type IntoIter = Iter<'a>;

    fn into_iter(self) -> Self::IntoIter {
        Iter {
            inner: self.hpo_terms.iter(),
            ontology: self,
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn add_terms_from_bytes() {
        let test_terms = [
            ("t1", 1u32),
            ("Term with a very long name", 2u32),
            ("", 3u32),
            ("Abnormality", 4u32),
        ];

        let mut ont = Ontology::default();

        let mut v: Vec<u8> = Vec::new();
        for (name, id) in test_terms {
            let t = HpoTermInternal::new(String::from(name), id.into());
            v.append(&mut t.as_bytes());
        }
        ont.add_terms_from_bytes(Bytes::new(&v, parser::binary::BinaryVersion::V1));
        assert_eq!(ont.len(), 4);
    }

    #[test]
    fn add_parents_from_bytes() {
        let test_terms = [
            ("t1", 1u32),
            ("Term with a very long name", 2u32),
            ("", 3u32),
            ("Abnormality", 4u32),
        ];

        let mut ont = Ontology::default();

        for (name, id) in test_terms {
            ont.add_term(HpoTermInternal::new(String::from(name), id.into()));
        }
        assert_eq!(ont.len(), 4);

        // The fake term has the same HpoTermId as one of of the Test ontology
        let mut fake_term = HpoTermInternal::new(String::new(), 3u32.into());
        fake_term.add_parent(1u32);
        fake_term.add_parent(2u32);

        let bytes = fake_term.parents_as_byte();

        ont.add_parent_from_bytes(&bytes[..]);

        assert_eq!(ont.get_unchecked(3u32).parents().len(), 2);
        assert_eq!(ont.get_unchecked(1u32).children().len(), 1);
        assert_eq!(ont.get_unchecked(2u32).children().len(), 1);
    }

    #[test]
    fn parse_hpo_version() {
        let mut ont = Ontology::default();

        // 7*256 + 231 == 2023
        let v = [7u8, 231u8, 1u8, 31u8];
        ont.hpo_version_from_bytes(&Bytes::new(&v, BinaryVersion::V2))
            .unwrap();
        assert_eq!(ont.hpo_version, (2023, 1, 31));
        assert_eq!(ont.hpo_version(), "2023-01-31");

        ont.hpo_version_from_bytes(&Bytes::new(&v, BinaryVersion::V1))
            .unwrap();
        assert_eq!(ont.hpo_version(), "0000-00-00");
    }

    #[test]
    fn check_v2_parsing() {
        let ont = Ontology::from_binary("tests/example.hpo").unwrap();

        assert_eq!(ont.hpo_version, (2023, 1, 31));
        assert_eq!(ont.hpo_version(), "2023-01-31");
    }

    #[test]
    fn compare_v1_v2() {
        let ont1 = Ontology::from_binary("tests/example_v1.hpo").unwrap();
        let ont2 = Ontology::from_binary("tests/example.hpo").unwrap();

        let diff = ont1.compare(&ont2);
        assert_eq!(diff.added_hpo_terms().len(), 0);
        assert_eq!(diff.removed_hpo_terms().len(), 0);
        assert_eq!(diff.changed_hpo_terms().len(), 0);

        assert_eq!(diff.added_genes().len(), 0);
        assert_eq!(diff.removed_genes().len(), 0);
        assert_eq!(diff.changed_genes().len(), 0);

        assert_eq!(diff.added_omim_diseases().len(), 0);
        assert_eq!(diff.removed_omim_diseases().len(), 0);
        assert_eq!(diff.changed_omim_diseases().len(), 0);
    }
}
