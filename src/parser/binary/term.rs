//! Create [`HpoTermInternal`] from binary representation
//! and vice versa

use super::Bytes;
use crate::term::internal::HpoTermInternal;
use crate::HpoError;

/// Creates an `HpoTermInternal` from bytes
///
/// The binary layout is defined as:
///
/// | Byte offset | Number of bytes | Description |
/// | --- | --- | --- |
/// | 0 | 4 | The total length of the binary data blob as big-endian `u32` |
/// | 4 | 4 | The Term ID as big-endian `u32` |
/// | 8 | 1 | The length of the Term Name (converted to a u8 vector) as a `u8` |
/// | 9 | n | The Term name as u8 vector. If the name has more than 255 bytes, it is trimmed to 255 |
///
/// # Panics
///
/// This method will panic if the total byte length is longer than `u32::MAX`
pub(crate) fn from_bytes_v1(bytes: Bytes) -> Result<HpoTermInternal, HpoError> {
    if bytes.len() < 4 + 4 + 1 {
        return Err(HpoError::ParseBinaryError);
    }
    let total_len = u32::from_be_bytes([bytes[0], bytes[1], bytes[2], bytes[3]]);

    let id = u32::from_be_bytes([bytes[4], bytes[5], bytes[6], bytes[7]]);
    let name_len = bytes[8] as usize;

    if bytes.len() < 4 + 4 + 1 + name_len {
        return Err(HpoError::ParseBinaryError);
    }

    let Ok(name) = String::from_utf8(bytes[9..total_len as usize].to_vec()) else {
        return Err(HpoError::ParseBinaryError);
    };
    Ok(HpoTermInternal::new(name, id.into()))
}

/// Creates an `HpoTermInternal` from bytes
///
/// The binary layout is defined as:
///
/// | Byte offset | Number of bytes | Description |
/// | --- | --- | --- |
/// | 0 | 4 | The total length of the binary data blob as big-endian `u32` |
/// | 4 | 4 | The Term ID as big-endian `u32` |
/// | 8 | 1 | The length of the Term Name (converted to a u8 vector) as a `u8` |
/// | 9 | n | The Term name as u8 vector. If the name has more than 255 bytes, it is trimmed to 255 |
/// | 9 + n | 1 | Flags, currently only obsolete
/// | 10 + n | 4 | Term ID of a replacement term as big-endian `u32` or 0 if None |
///
/// # Panics
///
/// This method will panic if the total byte length is longer than `u32::MAX`
pub(crate) fn from_bytes_v2(bytes: Bytes) -> Result<HpoTermInternal, HpoError> {
    if bytes.len() < 14 {
        return Err(HpoError::ParseBinaryError);
    }

    let id = u32::from_be_bytes([bytes[4], bytes[5], bytes[6], bytes[7]]);
    let name_len = bytes[8] as usize;

    if bytes.len() < 14 + name_len {
        return Err(HpoError::ParseBinaryError);
    }

    let Ok(name) = String::from_utf8(bytes[9..9 + name_len].to_vec()) else {
        return Err(HpoError::ParseBinaryError);
    };
    let mut term = HpoTermInternal::new(name, id.into());

    if bytes[9 + name_len] & 1u8 == 1u8 {
        *term.obsolete_mut() = true;
    }
    let offset = 10 + name_len;
    let replacement_id = u32::from_be_bytes([
        bytes[offset],
        bytes[offset + 1],
        bytes[offset + 2],
        bytes[offset + 3],
    ]);
    if replacement_id != 0 {
        *term.replacement_mut() = Some(replacement_id.into());
    }
    Ok(term)
}
