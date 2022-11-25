use crate::annotations::{Genes, GeneId};
use crate::term::{HpoChildren, HpoParents, HpoTermId};
use crate::HashSet;
use crate::DEFAULT_NUM_ALL_PARENTS;
use crate::DEFAULT_NUM_PARENTS;
use crate::DEFAULT_NUM_GENES;
use crate::term::InformationContent;

#[derive(Debug)]
pub struct HpoTermInternal {
    id: HpoTermId,
    name: String,
    parents: HpoParents,
    all_parents: HpoParents,
    children: HpoChildren,
    genes: Genes,
    ic: InformationContent,
}

impl HpoTermInternal {
    pub fn new(name: &str) -> HpoTermInternal {
        HpoTermInternal {
            id: name.into(),
            name: name.to_string(),
            parents: HashSet::with_capacity(DEFAULT_NUM_PARENTS),
            all_parents: HashSet::with_capacity(DEFAULT_NUM_ALL_PARENTS),
            children: HashSet::with_capacity(DEFAULT_NUM_PARENTS),
            genes: HashSet::with_capacity(DEFAULT_NUM_GENES),
            ic: InformationContent::default()
        }
    }

    pub fn id(&self) -> &HpoTermId {
        &self.id
    }

    pub fn name(&self) -> &str {
        &self.name
    }

    pub fn parents(&self) -> &HpoParents {
        &self.parents
    }

    pub fn all_parents(&self) -> &HpoParents {
        &self.all_parents
    }

    pub fn all_parents_mut(&mut self) -> &mut HpoParents {
        &mut self.all_parents
    }

    pub fn genes(&self) -> &Genes {
        &self.genes
    }

    pub fn parents_cached(&self) -> bool {
        if self.parents.is_empty() {
            true
        } else {
            !self.all_parents.is_empty()
        }
    }

    pub fn add_parent(&mut self, parent_id: HpoTermId) {
        self.parents.insert(parent_id);
    }

    pub fn add_child(&mut self, child_id: HpoTermId) {
        self.children.insert(child_id);
    }

    pub fn add_gene(&mut self, gene_id: GeneId) -> bool {
        self.genes.insert(gene_id)
    }
}


impl PartialEq for HpoTermInternal {
    fn eq(&self, other: &Self) -> bool {
        self.id == other.id
    }
}

impl Eq for HpoTermInternal {}
