#![deny(clippy::pedantic)]
#![allow(clippy::must_use_candidate)]
#![allow(clippy::module_name_repetitions)]
// #![warn(missing_docs)]
#![warn(missing_doc_code_examples)]
#![doc = include_str!("../README.md")]

use core::fmt::Debug;
use std::num::ParseIntError;
use thiserror::Error;

pub mod annotations;
mod matrix;
mod ontology;
mod parser;
mod set;
pub mod similarity;
pub mod stats;
pub mod term;

pub use ontology::Ontology;
pub use set::HpoSet;
pub use term::{HpoTerm, HpoTermId};

const DEFAULT_NUM_PARENTS: usize = 10;
const DEFAULT_NUM_ALL_PARENTS: usize = 50;
const DEFAULT_NUM_GENES: usize = 50;
const DEFAULT_NUM_OMIM: usize = 20;
const MAX_HPO_ID_INTEGER: usize = 10_000_000;

const OBO_FILENAME: &str = "hp.obo";
const GENE_FILENAME: &str = "phenotype_to_genes.txt";
const DISEASE_FILENAME: &str = "phenotype.hpoa";

#[derive(Error, Debug)]
/// Main Error type for this crate
pub enum HpoError {
    /// Indicates that a method or feature is not yet implemented
    #[error("not implemented")]
    NotImplemented,
    /// The term does not exist in the Ontology
    #[error("term does not exist")]
    DoesNotExist,
    /// Parsing of an integer failed
    #[error("unable to parse Integer")]
    ParseIntError,
    /// An error occured during parsing the binary HPO data
    #[error("unable to parse binary data")]
    ParseBinaryError,
    /// Opening the file was not able - check if the file is present
    #[error("cannot open file {0}")]
    CannotOpenFile(String),
    /// Failed to convert an integer to a float
    #[error("cannot convert int to float")]
    TryFromIntError(#[from] std::num::TryFromIntError),
}

impl From<ParseIntError> for HpoError {
    fn from(_: ParseIntError) -> Self {
        HpoError::ParseIntError
    }
}

/// Shortcut for `Result<T, HpoError>`
pub type HpoResult<T> = Result<T, HpoError>;

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

fn f32_from_usize(n: usize) -> HpoResult<f32> {
    let intermediate: u16 = n.try_into()?;
    Ok(intermediate.into())
}
