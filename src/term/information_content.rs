#[derive(Default, Debug)]
pub struct InformationContent {
    gene: f32,
    omim: f32,
}

impl InformationContent {
    pub fn gene(&self) -> f32 {
        self.gene
    }

    pub fn gene_mut(&mut self) -> &mut f32 {
        &mut self.gene
    }

    pub fn omim_disease(&self) -> f32 {
        self.omim
    }

    pub fn omim_disease_mut(&mut self) -> &mut f32 {
        &mut self.omim
    }

    pub fn get_kind(&self, kind: &InformationContentKind) -> f32 {
        match kind {
            InformationContentKind::Gene => self.gene(),
            InformationContentKind::Omim => self.omim_disease(),
        }
    }
}

pub enum InformationContentKind {
    Gene,
    Omim,
}
