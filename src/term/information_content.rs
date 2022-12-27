/// The quality (in terms of specificity) of an HPO term
///
/// Information content describes how specific an HPO term is
/// in relation to its associated genes or diseases
///
/// For more information, see
/// Resnik P: Using information content to evaluate semantic similarity in a taxonomy.
/// Proceedings of the 14th International Joint Conference on Artificial Intelligence:
/// August 20-25. 1995, Morgan Kaufmann, San Francisco CA: Montreal, Canada
#[derive(Default, Debug)]
pub struct InformationContent {
    gene: f32,
    omim: f32,
}

impl InformationContent {
    /// The Gene-specific information content
    pub fn gene(&self) -> f32 {
        self.gene
    }

    /// A mutable reference to the Gene-specific information content
    pub fn gene_mut(&mut self) -> &mut f32 {
        &mut self.gene
    }

    /// The OMIM-disease-specific information content
    pub fn omim_disease(&self) -> f32 {
        self.omim
    }

    /// A mutable reference to the OMIM-disease-specific information content
    pub fn omim_disease_mut(&mut self) -> &mut f32 {
        &mut self.omim
    }

    /// Returns the information content of the provided kind
    pub fn get_kind(&self, kind: &InformationContentKind) -> f32 {
        match kind {
            InformationContentKind::Gene => self.gene(),
            InformationContentKind::Omim => self.omim_disease(),
        }
    }
}

/// Different types of information contents
pub enum InformationContentKind {
    /// Information content related to the associated genes
    Gene,
    /// Information content related to the associated OMIM-diseases
    Omim,
}
