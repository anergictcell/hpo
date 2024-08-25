//! Parsing the ontology
use crate::parser::binary::Bytes;
use crate::{HpoError, HpoResult};

/// Identifies the encoding version of the provided binary data
///
/// # Errors
///
/// This method can fail for the following reasons:
/// - Incorrect byte encoding: [`HpoError::ParseBinaryError`]
/// - Invalid version: [`HpoError::NotImplemented`]
///
pub(crate) fn version(bytes: &[u8]) -> HpoResult<Bytes> {
    if bytes.len() < 5 {
        return Err(HpoError::ParseBinaryError);
    }

    if bytes[0..3] == [0x48, 0x50, 0x4f] {
        match bytes[3] {
            3u8 => Ok(Bytes::new(&bytes[4..], super::BinaryVersion::V3)),
            2u8 => Ok(Bytes::new(&bytes[4..], super::BinaryVersion::V2)),
            _ => Err(HpoError::NotImplemented),
        }
    } else {
        Ok(Bytes::new(bytes, super::BinaryVersion::V1))
    }
}
