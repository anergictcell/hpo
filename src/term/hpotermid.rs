use core::fmt::Debug;
use std::fmt::Display;

use crate::{OntologyResult, HpoError};

#[derive(Copy, Clone, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub struct HpoTermId {
    inner: u32,
}

impl HpoTermId {
    fn new(s: &str) -> HpoTermId {
        HpoTermId {
            inner: s[3..].parse::<u32>().expect("Invalid HPO-Term ID"),
        }
    }

    pub fn to_usize(&self) -> usize {
        self.inner.try_into().unwrap()
    }
}

impl TryFrom<&str> for HpoTermId {
    type Error = HpoError;
    fn try_from(s: &str) -> OntologyResult<Self> {
        Ok(HpoTermId {inner: s[3..].parse::<u32>()?})
    }
}

// impl From<&str> for HpoTermId {
//     fn from(s: &str) -> Self {
//         HpoTermId::new(s)
//     }
// }

impl From<String> for HpoTermId {
    fn from(s: String) -> Self {
        HpoTermId::new(&s)
    }
}

impl From<usize> for HpoTermId {
    fn from(n: usize) -> Self {
        Self {
            inner: n.try_into().unwrap(),
        }
    }
}

impl Debug for HpoTermId {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "HpoTermId({})", self)
    }
}

impl Display for HpoTermId {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "HP:{:07}", self.inner)
    }
}

impl PartialEq<str> for HpoTermId {
    fn eq(&self, other: &str) -> bool {
        self == &HpoTermId::new(other)
    }
}
