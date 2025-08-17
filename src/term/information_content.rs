use crate::{f32_from_usize, HpoResult};

/// The quality (in terms of specificity) of an HPO term
///
/// Information content describes how specific an HPO term is
/// in relation to its associated genes or diseases
///
/// For more information, see
/// Resnik P: Using information content to evaluate semantic similarity in a taxonomy.
/// Proceedings of the 14th International Joint Conference on Artificial Intelligence:
/// August 20-25. 1995, Morgan Kaufmann, San Francisco CA: Montreal, Canada
#[derive(Default, Debug, Clone)]
pub struct InformationContent {
    gene: f32,
    omim: f32,
    orpha: f32,
    custom: f32,
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

    /// The ORPHA-disease-specific information content
    pub fn orpha_disease(&self) -> f32 {
        self.orpha
    }

    /// A mutable reference to the ORPHA-disease-specific information content
    pub fn orpha_disease_mut(&mut self) -> &mut f32 {
        &mut self.orpha
    }

    /// Returns the information content of the provided kind
    pub fn get_kind(&self, kind: &InformationContentKind) -> f32 {
        match kind {
            InformationContentKind::Gene => self.gene(),
            InformationContentKind::Omim => self.omim_disease(),
            InformationContentKind::Orpha => self.orpha_disease(),
            InformationContentKind::Custom => self.custom(),
        }
    }

    fn calculate(total: usize, current: usize) -> HpoResult<f32> {
        if total == 0 || current == 0 {
            return Ok(0.0);
        }
        let total = f32_from_usize(total)?;
        let current = f32_from_usize(current)?;

        Ok((current / total).ln() * -1.0)
    }

    /// Calculates and caches the gene `InformationContent`
    ///
    /// # Errors
    ///
    /// This method returns an error if there are more genes than `u16::MAX`
    /// because larger numbers can't be safely converted to `f32`
    pub fn set_gene(&mut self, total: usize, current: usize) -> HpoResult<()> {
        self.gene = Self::calculate(total, current)?;
        Ok(())
    }

    /// Calculates and caches the OMIM `InformationContent`
    ///
    /// # Errors
    ///
    /// This method returns an error if there are more genes than `u16::MAX`
    /// because larger numbers can't be safely converted to `f32`
    pub fn set_omim_disease(&mut self, total: usize, current: usize) -> HpoResult<()> {
        self.omim = Self::calculate(total, current)?;
        Ok(())
    }

    /// Calculates and caches the ORPHA `InformationContent`
    ///
    /// # Errors
    ///
    /// This method returns an error if there are more ORPHA diseases than `u16::MAX`
    /// because larger numbers can't be safely converted to `f32`
    pub fn set_orpha_disease(&mut self, total: usize, current: usize) -> HpoResult<()> {
        self.orpha = Self::calculate(total, current)?;
        Ok(())
    }

    /// A custom Information content that can be provided by a client
    /// based on their own method
    pub fn custom(&self) -> f32 {
        self.custom
    }

    /// A mutable reference to the custom information content
    pub fn custom_mut(&mut self) -> &mut f32 {
        &mut self.custom
    }
}

/// Different types of information contents
#[derive(Debug, Copy, Clone, PartialEq, Eq, Hash)]
pub enum InformationContentKind {
    /// Information content related to the associated genes
    Gene,
    /// Information content related to the associated OMIM-diseases
    Omim,
    /// Information content related to the associated ORPHA-diseases
    Orpha,
    /// This information content can be custom defined by clients
    Custom,
}
