use std::collections::HashSet;


use crate::HpoTermId;
use crate::HpoError;
use crate::Ontology;
use crate::term::HpoGroup;



pub type OmimDiseases = HashSet<OmimDiseaseId>;

#[derive(Clone, Copy, Default, Debug, Hash, PartialEq, Eq)]
pub struct OmimDiseaseId {
    inner: usize
}

impl TryFrom<&str> for OmimDiseaseId {
    type Error = HpoError;
    fn try_from(value: &str) -> Result<Self, Self::Error> {
        Ok(OmimDiseaseId {inner: value.parse::<usize>()?})
    }
}


pub struct OmimDisease {
    id: OmimDiseaseId,
    name: String,
    hpos: HpoGroup
}

impl OmimDisease {
    pub fn new(id: OmimDiseaseId, name: &str) -> OmimDisease {
        Self { 
            name: name.to_string(),
            id,
            hpos: HpoGroup::default()
        }
    }

    pub fn id(&self) -> &OmimDiseaseId {
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

impl PartialEq for OmimDisease {
    fn eq(&self, other: &OmimDisease) -> bool {
        self.id == other.id
    }
}
impl Eq for OmimDisease {}


pub struct OmimDiseaseIterator<'a> {
    ontology: &'a Ontology,
    diseases: std::collections::hash_set::Iter<'a, OmimDiseaseId>,
}

impl<'a> OmimDiseaseIterator<'a> {
    pub fn new(diseases: &'a OmimDiseases, ontology: &'a Ontology) -> Self {
        OmimDiseaseIterator {
            diseases: diseases.iter(),
            ontology,
        }
    }
}

impl<'a> std::iter::Iterator for OmimDiseaseIterator<'a> {
    type Item = &'a OmimDisease;
    fn next(&mut self) -> Option<Self::Item> {
        match self.diseases.next() {
            Some(omim_id) => {
                Some(self.ontology.get_omim_disease(omim_id).unwrap())
            }
            None => None,
        }
    }
}
