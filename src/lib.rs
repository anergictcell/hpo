#![doc = include_str!("../README.md")]
#![allow(dead_code)]
use core::fmt::Debug;
use std::num::ParseIntError;
use thiserror::Error;

pub mod annotations;
mod matrix;
mod ontology;
pub mod parser;
mod set;
pub mod similarity;
pub mod term;

pub use ontology::Ontology;
pub use similarity::{GraphIc, Similarity};
pub use term::{HpoParents, HpoTerm, HpoTermId, InformationContentKind};

const DEFAULT_NUM_PARENTS: usize = 10;
const DEFAULT_NUM_ALL_PARENTS: usize = 50;
const DEFAULT_NUM_GENES: usize = 50;
const DEFAULT_NUM_OMIM: usize = 20;
const MAX_HPO_ID_INTEGER: usize = 10_000_000;

#[derive(Error, Debug)]
pub enum HpoError {
    #[error("not implemented")]
    NotImplemented,
    #[error("term does not exist")]
    DoesNotExist,
    #[error("unable to parse Integer")]
    ParseIntError,
}

impl From<ParseIntError> for HpoError {
    fn from(_: ParseIntError) -> Self {
        HpoError::ParseIntError
    }
}

type OntologyResult<T> = Result<T, HpoError>;
