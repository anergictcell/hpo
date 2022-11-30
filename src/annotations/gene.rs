use core::fmt::Debug;
use std::cmp::PartialEq;
use std::collections::HashSet;
use std::convert::TryFrom;

use crate::term::HpoGroup;
use crate::HpoError;
use crate::HpoTermId;
use crate::Ontology;
use crate::OntologyResult;

pub type Genes = HashSet<GeneId>;

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

#[derive(Default, Debug)]
pub struct Gene {
    id: GeneId,
    name: String,
    hpos: HpoGroup,
}

impl Gene {
    pub fn new(id: GeneId, name: &str) -> Gene {
        Gene {
            id,
            name: name.to_string(),
            hpos: HpoGroup::default(),
        }
    }
    pub fn from_parts(id: &str, name: &str) -> OntologyResult<Gene> {
        Ok(Gene {
            id: GeneId::try_from(id)?,
            name: name.to_string(),
            hpos: HpoGroup::default(),
        })
    }

    pub fn id(&self) -> &GeneId {
        &self.id
    }

    pub fn name(&self) -> &str {
        &self.name
    }

    pub fn symbol(&self) -> &str {
        &self.name
    }

    pub fn add_term(&mut self, term_id: HpoTermId) -> bool {
        self.hpos.insert(term_id)
    }
}

impl PartialEq for Gene {
    fn eq(&self, other: &Gene) -> bool {
        self.id == other.id
    }
}
impl Eq for Gene {}

pub struct GeneIterator<'a> {
    ontology: &'a Ontology,
    genes: std::collections::hash_set::Iter<'a, GeneId>,
}

impl<'a> GeneIterator<'a> {
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
            Some(gene_id) => Some(self.ontology.get_gene(gene_id).unwrap()),
            None => None,
        }
    }
}

impl Debug for GeneIterator<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "GeneIterator")
    }
}
