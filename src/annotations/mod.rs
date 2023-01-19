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
pub use gene::{Gene, GeneId, GeneIterator, Genes};

mod omim_disease;
pub use omim_disease::{OmimDisease, OmimDiseaseId, OmimDiseaseIterator, OmimDiseases};
