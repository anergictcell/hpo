#![allow(dead_code)]
use core::fmt::Debug;
use std::num::ParseIntError;
use thiserror::Error;

pub mod parser;
pub mod term;
pub mod annotations;
mod ontology;

pub use term::{HpoTerm, HpoTermId, HpoParents, InformationContentKind};
pub use ontology::Ontology;

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

enum HpoLink {
    Parent,
}


pub trait Similarity {
    fn calculate(&self, a: &HpoTerm, b: &HpoTerm) -> f32;
}

pub struct GraphIc {
    method: InformationContentKind
}

impl GraphIc {
    pub fn new(method: InformationContentKind) -> Self {
        Self {method}
    }
}

impl Similarity for GraphIc {
    fn calculate(&self, a: &HpoTerm, b: &HpoTerm) -> f32 {
        if a.id() == b.id() {
            return 1.0
        }

        let ic_union: f32 = a.union_ancestors(b)
            .map(|p| p.information_content().get_kind(&self.method))
            .sum();

        if ic_union == 0.0 {
            return 0.0
        }

        let ic_common: f32 = a.common_ancestors(b)
            .map(|p| p.information_content().get_kind(&self.method))
            .sum();


        ic_common/ic_union
    }
}