//! Genes and Diseases are linked to HPO terms and make up secondary annotations
//!
//! This module contains structs to represent [`Gene`]s and [`OmimDisease`]s
//! and iterators.
//!
//! The underlying principle for all annotations is the same:
//! - Each record (gene or disease) has a unique numerical identifier.
//! - Each record holds information about which HPO terms it is related to
//! - Records can be grouped in Sets
//!
//! For more information check out the individual sections for [`Gene`] and [`OmimDisease`]

mod gene;
use core::fmt::Debug;
use core::hash::Hash;
pub use gene::{Gene, GeneId, GeneIterator, Genes};
use std::fmt::Display;

mod disease;
mod omim_disease;
mod orpha_disease;
pub use disease::Disease;
pub use omim_disease::{
    OmimDisease, OmimDiseaseFilter, OmimDiseaseId, OmimDiseaseIterator, OmimDiseases,
};
pub use orpha_disease::{OrphaDisease, OrphaDiseaseId, OrphaDiseaseIterator, OrphaDiseases};

/// All annotations ([`Gene`]s, [`OmimDisease`]s) are defined by a unique ID that is constrained by this trait
///
/// The ID must be unique only within the annotation type, i.e. a gene and a disease
/// can have the same ID.
pub trait AnnotationId:
    Clone
    + Copy
    + Debug
    + Hash
    + PartialEq
    + PartialOrd
    + Eq
    + Ord
    + Display
    + From<u32>
    + for<'a> TryFrom<&'a str>
{
    /// Return the integer representation of the annotation ID
    fn as_u32(&self) -> u32;

    /// Returns the memory representation of the inner integer as a byte array in big-endian (network) byte order.
    fn to_be_bytes(&self) -> [u8; 4] {
        self.as_u32().to_be_bytes()
    }
}
