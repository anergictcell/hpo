use crate::annotations::{GeneId, Genes};
use crate::annotations::{OmimDiseaseId, OmimDiseases};
use crate::term::HpoGroup;
use crate::term::InformationContent;
use crate::term::{HpoChildren, HpoParents, HpoTermId};
use crate::DEFAULT_NUM_ALL_PARENTS;
use crate::DEFAULT_NUM_GENES;
use crate::DEFAULT_NUM_OMIM;
use crate::DEFAULT_NUM_PARENTS;

#[derive(Debug)]
pub struct HpoTermInternal {
    id: HpoTermId,
    name: String,
    parents: HpoParents,
    all_parents: HpoParents,
    children: HpoChildren,
    genes: Genes,
    omim_diseases: OmimDiseases,
    ic: InformationContent,
    obsolete: bool,
    replacement: Option<HpoTermId>,
}

impl Default for HpoTermInternal {
    fn default() -> Self {
        HpoTermInternal::new("HP:0000000")
    }
}

impl HpoTermInternal {
    pub fn new(name: &str) -> HpoTermInternal {
        HpoTermInternal {
            id: name.into(),
            name: name.to_string(),
            parents: HpoGroup::with_capacity(DEFAULT_NUM_PARENTS),
            all_parents: HpoGroup::with_capacity(DEFAULT_NUM_ALL_PARENTS),
            children: HpoChildren::with_capacity(DEFAULT_NUM_PARENTS),
            genes: Genes::with_capacity(DEFAULT_NUM_GENES),
            omim_diseases: OmimDiseases::with_capacity(DEFAULT_NUM_OMIM),
            ic: InformationContent::default(),
            obsolete: false,
            replacement: None,
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

    pub fn omim_diseases(&self) -> &OmimDiseases {
        &self.omim_diseases
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

    pub fn add_omim_disease(&mut self, omim_disease_id: OmimDiseaseId) -> bool {
        self.omim_diseases.insert(omim_disease_id)
    }

    pub fn information_content(&self) -> &InformationContent {
        &self.ic
    }

    pub fn information_content_mut(&mut self) -> &mut InformationContent {
        &mut self.ic
    }
}

impl PartialEq for HpoTermInternal {
    fn eq(&self, other: &Self) -> bool {
        self.id == other.id
    }
}

impl Eq for HpoTermInternal {}
