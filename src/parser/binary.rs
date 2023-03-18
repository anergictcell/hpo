//! The binary parser handles the parsing and generation of the Ontology in
//! binary format.
//! The idea is to always support old binary formats as input but only
//! generate the newest binary format.
//! That way we maintain backwards compatibility and allow clients
//! to update their versions to the newest one.
pub(crate) mod ontology;
pub(crate) mod term;
use std::fmt::Display;

use crate::{term::internal::HpoTermInternal, HpoError};

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub(crate) enum BinaryVersion {
    V1,
    V2,
}

impl TryFrom<u8> for BinaryVersion {
    type Error = HpoError;
    fn try_from(value: u8) -> Result<Self, Self::Error> {
        match value {
            1u8 => Ok(BinaryVersion::V1),
            2u8 => Ok(BinaryVersion::V2),
            _ => Err(HpoError::NotImplemented),
        }
    }
}

impl Display for BinaryVersion {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                BinaryVersion::V1 => "1",
                BinaryVersion::V2 => "2",
            }
        )
    }
}

#[derive(Debug, Clone, Copy)]
pub(crate) struct Bytes<'a> {
    data: &'a [u8],
    version: BinaryVersion,
}

impl<'a> Bytes<'a> {
    pub fn new(data: &'a [u8], version: BinaryVersion) -> Self {
        Self { data, version }
    }

    pub fn len(&self) -> usize {
        self.data.len()
    }

    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    pub fn subset<T>(&self, idx: T) -> Bytes<'_>
    where
        T: std::slice::SliceIndex<[u8], Output = [u8]>,
    {
        Bytes::new(&self.data[idx], self.version)
    }

    pub fn version(&self) -> BinaryVersion {
        self.version
    }

    pub fn u32_prefix(&self) -> u32 {
        assert!(self.len() > 4, "4 bytes are required to extract u32");
        u32::from_be_bytes([self[0], self[1], self[2], self[3]])
    }
}

impl<Idx> std::ops::Index<Idx> for Bytes<'_>
where
    Idx: std::slice::SliceIndex<[u8]>,
{
    type Output = Idx::Output;

    fn index(&self, idx: Idx) -> &Self::Output {
        &self.data[idx]
    }
}

/// Builder to crate multiple [`HpoTermInternal`] from raw bytes
pub(crate) struct BinaryTermBuilder<'a>(Bytes<'a>, usize);

impl<'a> BinaryTermBuilder<'a> {
    /// Crates a new [`BinaryTermBuilder`]
    pub fn new(bytes: Bytes<'a>) -> Self {
        Self(bytes, 0)
    }
}

impl Iterator for BinaryTermBuilder<'_> {
    type Item = HpoTermInternal;
    fn next(&mut self) -> Option<Self::Item> {
        let bytes = self.0.subset(self.1..);

        if bytes.is_empty() {
            return None;
        }

        let term_len = bytes.u32_prefix() as usize;

        assert!(
            (bytes.len() >= term_len),
            "Invalid bytes left over in BinaryTermBuilder"
        );

        self.1 += term_len;
        Some(HpoTermInternal::try_from(bytes).expect("Invalid byte input"))
    }
}
