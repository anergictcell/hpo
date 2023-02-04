//! Termarena is a private module that contains the raw data model holding
//! all HPO-Terms.
//!
//! The data layout can be described as having two `Vec`:
//!
//! - `terms` holds all `HpoTerm`s (`HpoTermInternal`) with one element for every term
//!   in the ontology. The terms are not sorted and are added in the same order they
//!   are defined in the source data. The order has no impact.
//! - `ids` is a lookup table that contains one element for every possible [`HpoTermId`].
//!   HPO-Term-IDs can be any integer between 1 and 10 Million, so it has 10,000,000 elements.
//!   Each element contains either `0` (no such term present) or the index of the
//!   corresponding term in `terms`.
//!
//! To account for the fact that a `0` value in `ids` signals that no term is present, the
//! first entry in `terms` is a fake `HpoTerm`.
//!
//! # Data layout example
//!
//! ## `ids` Vector
//!
//! |  i  | value |
//! | --- | --- |
//! |  0  |  0  |
//! |  1  |  1  |
//! |  2  |  0  |
//! |  3  |  5  |
//! |  4  |  2  |
//! |  5  |  0  |
//! |  6  |  3  |
//! |  7  |  4  |
//!
//!
//! ## `terms` Vector
//!
//! |  i  | Term-ID |
//! | --- | --- |
//! |  0  |  0 (Fake)  |
//! |  1  |  1  |
//! |  2  |  4  |
//! |  3  |  6  |
//! |  4  |  7  |
//! |  5  |  3  |
//!
//!
//! This is a bit like a cheap version of a Hashmap. The advantage over a hashmap are
//!
//! - very cheap hash algorithm - The integer representation for each term is given with the term ID
//! - we don't have to worry about or account for hash collisions
//! - we know the size of the lookup vector and can allocate it once
//! - we have an extremely good estimate of the amount of HPO terms, so we can pre-allocate that vector as well
//! - every lookup is O(3)
//!     - get `idx` from `ids` vector
//!     - check if `idx` is 0
//!     - get `HpoTermInternal` from `terms` vector
//!
//! The disadvantages are
//!
//! - I built it myself and I don't have much prior experience with efficient hashmap data structures
//! - I'm sure better Hashmap implementations exist that would solve this problem in a more performant way
//!
//! I initially used the default `Hashmap` instead, but the performance with the custom lookup is significantly better.
//! Since the retrieval of `HpoTerm`s from the ontology is the most frequent operation, I wanted to optimize
//! this lookup as much as possible.
//!
//! If you're looking at this and are shaking your head, please let me know how to improve this. I'd love to hear feedback,
//! this part of the crate is the most performance-critical and I am also very interested in learning how to improve
//! such data structures.
//!

#![allow(clippy::slow_vector_initialization)]
use log::trace;
use log::warn;

use crate::term::internal::HpoTermInternal;
use crate::HpoTermId;

use crate::MAX_HPO_ID_INTEGER as HPO_TERM_NUMBERS;

pub(crate) struct Arena {
    terms: Vec<HpoTermInternal>,
    ids: Vec<usize>,
}

impl Default for Arena {
    fn default() -> Self {
        let mut ids = Vec::with_capacity(HPO_TERM_NUMBERS);
        ids.resize(HPO_TERM_NUMBERS, 0);
        let mut s = Self {
            terms: Vec::with_capacity(18_000),
            ids,
        };
        s.terms.push(HpoTermInternal::default());
        s
    }
}

impl Arena {
    /// Returns how many `HpoTerm`s are present in the Arena
    pub fn len(&self) -> usize {
        self.terms.len() - 1
    }

    /// Inserts a new [`HpoTerm`](crate::term::HpoTerm) into the arena
    ///
    /// If a term with the same ID exists already, it does nothing
    pub fn insert(&mut self, term: HpoTermInternal) {
        let id = term.id().to_usize();
        if self.ids[id] == 0 {
            let idx = self.terms.len();
            self.terms.push(term);
            self.ids[id] = idx;
        }
    }

    /// Returns the [`HpoTermInternal`] with the given `HpoTermId`
    ///
    /// If no such term is present, returns `None`
    pub fn get(&self, id: HpoTermId) -> Option<&HpoTermInternal> {
        match self.ids.get(id.to_usize()) {
            Some(0) => {
                trace!("Term does not exist in Arena: {}", id);
                None
            }
            Some(n) => Some(&self.terms[*n]),
            None => {
                warn!("Index of Arena out of bounds for {id}");
                None
            }
        }
    }

    /// Returns the [`HpoTermInternal`] with the given `HpoTermId`
    ///
    /// This method does not check if the term is actually present
    /// and should only ne used if you are sure that the term exists
    ///
    /// # Panics
    ///
    /// If no `HpoTerm` with the given ID exists in the Arena
    pub fn get_unchecked(&self, id: HpoTermId) -> &HpoTermInternal {
        &self.terms[self.ids[id.to_usize()]]
    }

    /// Returns a mutable reference to the [`HpoTermInternal`] with the given `HpoTermId`
    ///
    /// This method does not check if the term is actually present
    /// and should only ne used if you are sure that the term exists
    ///
    /// # Panics
    ///
    /// If no `HpoTerm` with the given ID exists in the Arena
    pub fn get_unchecked_mut(&mut self, id: HpoTermId) -> &mut HpoTermInternal {
        &mut self.terms[self.ids[id.to_usize()]]
    }

    /// Returns a mutable reference to the [`HpoTermInternal`] with the given `HpoTermId`
    ///
    /// If no such term is present, returns `None`
    pub fn get_mut(&mut self, id: HpoTermId) -> Option<&mut HpoTermInternal> {
        match self.ids.get(id.to_usize()) {
            Some(0) => None,
            Some(n) => Some(&mut self.terms[*n]),
            None => {
                warn!("Index of Arena out of bounds for {id}");
                None
            }
        }
    }

    /// Returns a slice of all [`HpoTermInternal`]
    ///
    /// This can be used for iterating all terms
    pub fn values(&self) -> &[HpoTermInternal] {
        &self.terms[1..]
    }

    /// Returns a mutable slice of all [`HpoTermInternal`]
    ///
    /// This can be used for iterating and modifying all terms
    pub fn values_mut(&mut self) -> &mut [HpoTermInternal] {
        &mut self.terms[1..]
    }

    /// Returns all [`HpoTermId`]s
    pub fn keys(&mut self) -> Vec<HpoTermId> {
        self.terms[1..].iter().map(|term| *term.id()).collect()
    }
}
