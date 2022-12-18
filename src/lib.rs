// #![warn(missing_docs)]
// #![warn(missing_doc_code_examples)]

#![doc = include_str!("../README.md")]
use core::fmt::Debug;
use std::num::ParseIntError;
use thiserror::Error;

pub mod annotations;
mod matrix;
mod ontology;
pub mod parser;
mod set;
pub mod similarity;
mod term;

pub use ontology::Ontology;
pub use similarity::{GraphIc, Similarity};
pub use term::{HpoParents, HpoTerm, HpoTermId, HpoTermIterator};
pub use term::{InformationContent, InformationContentKind};

const DEFAULT_NUM_PARENTS: usize = 10;
const DEFAULT_NUM_ALL_PARENTS: usize = 50;
const DEFAULT_NUM_GENES: usize = 50;
const DEFAULT_NUM_OMIM: usize = 20;
const MAX_HPO_ID_INTEGER: usize = 10_000_000;

const OBO_FILENAME: &str = "hp.obo";
const GENE_FILENAME: &str = "phenotype_to_genes.txt";
const DISEASE_FILENAME: &str = "phenotype.hpoa";

#[derive(Error, Debug)]
pub enum HpoError {
    #[error("not implemented")]
    NotImplemented,
    #[error("term does not exist")]
    DoesNotExist,
    #[error("unable to parse Integer")]
    ParseIntError,
    #[error("unable to parse binary data")]
    ParseBinaryError,
}

impl From<ParseIntError> for HpoError {
    fn from(_: ParseIntError) -> Self {
        HpoError::ParseIntError
    }
}

type OntologyResult<T> = Result<T, HpoError>;

/// Returns a u32 from the 4 first 4 bytes
///
/// This function is used in parsing binary ontology data
/// where u32 are used frequently. It is not recommended
/// to use this function elsewhere.
///
/// The u32 sould be big-endian encoded
fn u32_from_bytes(bytes: &[u8]) -> u32 {
    u32::from_be_bytes([bytes[0], bytes[1], bytes[2], bytes[3]])
}
