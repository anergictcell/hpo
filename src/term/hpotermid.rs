use core::fmt::Debug;
use std::fmt::Display;

use crate::{HpoError, HpoResult};

/// The ID of an HPO-Term (e.g. `HP:0000123`)
///
/// The ID must be convertible to an integer
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

    /// Convert `self` to a `usize`
    ///
    /// # Panics
    ///
    /// This method can panic on systems where `usize::MAX`
    /// is smaller than the `HpoTermId`
    pub fn to_usize(&self) -> usize {
        self.inner.try_into().unwrap()
    }

    /// Convert `self` to `u32`
    pub fn as_u32(&self) -> u32 {
        self.inner
    }

    /// Byte representation (big-endian) of self
    pub fn to_be_bytes(&self) -> [u8; 4] {
        self.inner.to_be_bytes()
    }
}

impl TryFrom<&str> for HpoTermId {
    type Error = HpoError;
    fn try_from(s: &str) -> HpoResult<Self> {
        Ok(HpoTermId {
            inner: s[3..].parse::<u32>()?,
        })
    }
}

impl From<String> for HpoTermId {
    fn from(s: String) -> Self {
        HpoTermId::new(&s)
    }
}

impl From<u16> for HpoTermId {
    fn from(n: u16) -> Self {
        Self { inner: n.into() }
    }
}

impl From<u32> for HpoTermId {
    fn from(inner: u32) -> Self {
        Self { inner }
    }
}

impl From<u64> for HpoTermId {
    fn from(n: u64) -> Self {
        Self {
            inner: n.try_into().unwrap(),
        }
    }
}

impl From<usize> for HpoTermId {
    fn from(n: usize) -> Self {
        Self {
            inner: n.try_into().unwrap(),
        }
    }
}

impl From<[u8; 4]> for HpoTermId {
    fn from(bytes: [u8; 4]) -> Self {
        Self {
            inner: u32::from_be_bytes(bytes),
        }
    }
}

impl Debug for HpoTermId {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "HpoTermId({self})")
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
