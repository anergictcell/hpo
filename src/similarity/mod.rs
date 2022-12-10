use crate::matrix::Matrix;
use crate::set::HpoSet;
use crate::HpoTerm;
use crate::InformationContentKind;

pub trait Similarity {
    fn calculate(&self, a: &HpoTerm, b: &HpoTerm) -> f32;
}

pub struct GraphIc {
    method: InformationContentKind,
}

impl GraphIc {
    pub fn new(method: InformationContentKind) -> Self {
        Self { method }
    }
}

impl Similarity for GraphIc {
    fn calculate(&self, a: &HpoTerm, b: &HpoTerm) -> f32 {
        if a.id() == b.id() {
            return 1.0;
        }

        let ic_union: f32 = a
            .union_ancestors(b)
            .map(|p| p.information_content().get_kind(&self.method))
            .sum();

        if ic_union == 0.0 {
            return 0.0;
        }

        let ic_common: f32 = a
            .common_ancestors(b)
            .map(|p| p.information_content().get_kind(&self.method))
            .sum();

        ic_common / ic_union
    }
}

pub trait SimilarityCombiner {
    fn combine(&self, m: &Matrix<f32>) -> f32;

    fn calculate(&self, m: &Matrix<f32>) -> f32 {
        if m.is_empty() {
            return 0.0;
        }
        self.combine(m)
    }

    fn row_maxes(&self, m: &Matrix<f32>) -> Vec<f32> {
        m.rows()
            .map(|row| {
                // I have no idea why, but I could not get a.max(b) to work
                // with the borrow checker
                row.reduce(|a, b| if a > b { a } else { b }).unwrap()
            })
            .cloned()
            .collect()
    }

    fn col_maxes(&self, m: &Matrix<f32>) -> Vec<f32> {
        m.cols()
            .map(|col| {
                // I have no idea why, but I could not get a.max(b) to work
                // with the borrow checker
                col.reduce(|a, b| if a > b { a } else { b }).unwrap()
            })
            .cloned()
            .collect()
    }

    fn dim_f32(&self, m: &Matrix<f32>) -> (f32, f32) {
        let (rows, cols) = m.dim();
        (usize_to_f32(rows), usize_to_f32(cols))
    }
}

pub enum StandardCombiner {
    FunSimAvg,
    FunSimMax,
    Bwa,
}

impl Default for StandardCombiner {
    fn default() -> Self {
        Self::FunSimAvg
    }
}

impl StandardCombiner {
    fn fun_sim_avg(&self, m: &Matrix<f32>) -> f32 {
        let (rows, cols) = self.dim_f32(m);
        let row_maxes = self.row_maxes(m);
        let col_maxes = self.col_maxes(m);
        let mut nom = row_maxes.iter().sum::<f32>() / rows;
        nom += col_maxes.iter().sum::<f32>() / cols;

        nom / 2.0
    }

    fn fun_sim_max(&self, m: &Matrix<f32>) -> f32 {
        let (rows, cols) = self.dim_f32(m);
        let row_maxes = self.row_maxes(m);
        let col_maxes = self.col_maxes(m);

        (row_maxes.iter().sum::<f32>() / rows).max(col_maxes.iter().sum::<f32>() / cols)
    }

    fn bwa(&self, m: &Matrix<f32>) -> f32 {
        let (rows, cols) = self.dim_f32(m);
        let row_maxes = self.row_maxes(m);
        let col_maxes = self.col_maxes(m);

        (row_maxes.iter().sum::<f32>() + col_maxes.iter().sum::<f32>()) / (rows + cols)
    }
}

impl SimilarityCombiner for StandardCombiner {
    fn combine(&self, m: &Matrix<f32>) -> f32 {
        match self {
            StandardCombiner::FunSimAvg => self.fun_sim_avg(m),
            StandardCombiner::FunSimMax => self.fun_sim_max(m),
            StandardCombiner::Bwa => self.bwa(m),
        }
    }
}

pub struct GroupSimilarity<T, C> {
    combiner: C,
    similarity: T,
}

impl<T: Similarity, C: SimilarityCombiner> GroupSimilarity<T, C> {
    pub fn new(combiner: C, similarity: T) -> Self {
        Self {
            combiner,
            similarity,
        }
    }

    fn calculate(&self, a: &HpoSet, b: &HpoSet) -> f32 {
        let mut v = Vec::with_capacity(a.len() * b.len());
        for t1 in a.into_iter() {
            for t2 in b.into_iter() {
                v.push(self.similarity.calculate(t1, t2))
            }
        }
        let m = Matrix::new(a.len(), b.len(), &v);
        self.combiner.calculate(&m)
    }
}

impl Default for GroupSimilarity<GraphIc, StandardCombiner> {
    fn default() -> Self {
        Self {
            combiner: StandardCombiner::default(),
            similarity: GraphIc::new(InformationContentKind::Omim),
        }
    }
}

/// This is a really weird way of converting a usize into a float but I
/// want to make sure the app crashes, so I don't want to use `as`.
fn usize_to_f32(n: usize) -> f32 {
    <usize as TryInto<u16>>::try_into(n)
        .expect("Matrix too large")
        .into()
}