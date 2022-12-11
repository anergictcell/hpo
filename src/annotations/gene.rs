use core::fmt::Debug;
use std::cmp::PartialEq;
use std::collections::HashSet;
use std::convert::TryFrom;
use std::fmt::Display;

use crate::term::HpoGroup;
use crate::HpoError;
use crate::HpoTermId;
use crate::Ontology;
use crate::OntologyResult;

/// A set of genes
///
/// The set does not contain [`Gene`]s itself, but only their [`GeneId`]s.
/// Currently implemented using [`HashSet`] but any other implementation
/// should work as well given that each GeneId must appear only once
/// and it provides an iterator of [`GeneId`]
pub type Genes = HashSet<GeneId>;


/// A unique identifier for a [`Gene`]
///
/// This value can - in theory - represent any numerical unique value.
/// When using the default JAX provided masterdata, it represents
/// the NCBI Gene ID
#[derive(Clone, Copy, Default, Debug, Hash, PartialEq, Eq)]
pub struct GeneId {
    inner: usize,
}

impl TryFrom<&str> for GeneId {
    type Error = HpoError;
    fn try_from(value: &str) -> Result<Self, Self::Error> {
        Ok(GeneId {
            inner: value.parse::<usize>()?,
        })
    }
}

impl Display for GeneId {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Gene NCBI:{}", self.inner)
    }
}


/// A single gene
///
/// A gene has a unique [`GeneId`] and a name (symbol) and is
/// connected to a set of HPO terms
#[derive(Default, Debug)]
pub struct Gene {
    id: GeneId,
    name: String,
    hpos: HpoGroup,
}

impl Gene {
    /// Initializes a new Gene
    ///
    /// This method should rarely, if ever, be used directly. The
    /// preferred way to create new genes is through [`Ontology::add_gene`]
    /// to ensure that each gene exists only once.
    pub fn new(id: GeneId, name: &str) -> Gene {
        Gene {
            id,
            name: name.to_string(),
            hpos: HpoGroup::default(),
        }
    }

    /// Initializes a new Gene from `str` values
    ///
    /// This method should rarely, if ever, be used directly. The
    /// preferred way to create new genes is through [`Ontology::add_gene`]
    /// to ensure that each gene exists only once.
    pub fn from_parts(id: &str, name: &str) -> OntologyResult<Gene> {
        Ok(Gene {
            id: GeneId::try_from(id)?,
            name: name.to_string(),
            hpos: HpoGroup::default(),
        })
    }

    /// The unique [`GeneId`] of the gene, most likely the NCBI Gene ID
    pub fn id(&self) -> &GeneId {
        &self.id
    }

    /// The name of the gene (gene symbol)
    pub fn name(&self) -> &str {
        &self.name
    }

    /// The gene symbol (identical to [`Gene::id`])
    pub fn symbol(&self) -> &str {
        &self.name
    }

    /// Connect another [HPO term](`HpoTerm`) to the gene 
    pub fn add_term(&mut self, term_id: HpoTermId) -> bool {
        self.hpos.insert(term_id)
    }

    /// The set of connected HPO terms
    pub fn hpo_terms(&self) -> &HpoGroup {
        &self.hpos
    }
}

impl PartialEq for Gene {
    fn eq(&self, other: &Gene) -> bool {
        self.id == other.id
    }
}
impl Eq for Gene {}

/// [`Gene`] Iterator 
pub struct GeneIterator<'a> {
    ontology: &'a Ontology,
    genes: std::collections::hash_set::Iter<'a, GeneId>,
}

impl<'a> GeneIterator<'a> {
    /// Initialize a new [`GeneIterator`]
    ///
    /// This method requires the [`Ontology`] as a parameter since
    /// the actual [`Gene`] entities are stored in it and not in [`Genes`]
    /// itself
    pub fn new(genes: &'a Genes, ontology: &'a Ontology) -> Self {
        GeneIterator {
            genes: genes.iter(),
            ontology,
        }
    }
}

impl<'a> std::iter::Iterator for GeneIterator<'a> {
    type Item = &'a Gene;
    fn next(&mut self) -> Option<Self::Item> {
        match self.genes.next() {
            Some(gene_id) => Some(self.ontology.gene(gene_id).unwrap()),
            None => None,
        }
    }
}

impl Debug for GeneIterator<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "GeneIterator")
    }
}
