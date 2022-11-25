#[derive(Default, Debug, )]
pub struct InformationContent {
    gene: f32,
    omim: f32
}

impl InformationContent {
    pub fn gene(&self) -> f32 {
        self.gene
    }

    pub fn gene_mut(&mut self) -> &mut f32 {
        &mut self.gene
    }

    pub fn omim(&self) -> f32 {
        self.omim
    }

    pub fn omim_mut(&mut self) -> &mut f32 {
        &mut self.omim
    }
}