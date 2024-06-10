use std::cmp::Ordering;
use std::collections::hash_map;
use std::collections::HashMap;

use crate::set::HpoSet;
use crate::utils::Combinations;
pub mod cluster;

use cluster::Cluster;
use cluster::ClusterVec;

#[derive(Debug, Default)]
struct DistanceMatrix(HashMap<(usize, usize), f32>);

impl<'a> DistanceMatrix {
    fn iter(&'a self) -> hash_map::Iter<(usize, usize), f32> {
        self.0.iter()
    }

    fn insert(&mut self, k: (usize, usize), v: f32) -> Option<f32> {
        self.0.insert(k, v)
    }

    fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    fn retain<F>(&mut self, f: F)
    where
        F: FnMut(&(usize, usize), &mut f32) -> bool,
    {
        self.0.retain(f);
    }

    fn get(&self, k: &(usize, usize)) -> Option<&f32> {
        self.0.get(k)
    }
}

/// Linkage matrices from `HpoSet`s
///
/// Crate a linkage matrix from a list of `HpoSet`s to use in dendograms
/// or other hierarchical cluster analyses
///
/// Provided algorithms for clustering
///
/// - [`Linkage::union`](`Linkage::union`): Create a new `HpoSet` for each cluster based on the union of
///   both combined clusters. This method becomes slow with growing input data
/// - [`Linkage::single`](`Linkage::single`): The minimum distance of each cluster's nodes to the other
///   nodes is used as distance for newly formed clusters. This is also known as the Nearest Point Algorithm.
/// - [`Linkage::complete`](`Linkage::complete`): The maximum distance of each cluster's nodes to the other
///   nodes is used as distance for newly formed clusters. This is also known by the Farthest Point Algorithm
///   or Voor Hees Algorithm.
/// - [`Linkage::average`](`Linkage::average`): The mean distance of each cluster's nodes to the other
///   nodes is used as distance for newly formed clusters. This is also called the UPGMA algorithm.
///
/// # Examples
///
/// ```rust
///
/// use hpo::Ontology;
/// use hpo::HpoSet;
/// use hpo::similarity::GroupSimilarity;
/// use hpo::utils::Combinations;
/// use hpo::stats::Linkage;
///
/// // This method can and should utilize parallel processing, e.g.
/// // using rayon iterators
/// fn distance(combs: Combinations<HpoSet<'_>>) -> Vec<f32> {
///     let sim = GroupSimilarity::default();
///     combs.map(|comp| {
///         1.0 - sim.calculate(comp.0, comp.1)
///     }).collect()
/// }
///
/// let ontology = Ontology::from_binary("tests/example.hpo").unwrap();
/// let sets = vec![
///     ontology.gene_by_name("KRAS").unwrap().to_hpo_set(&ontology),
///     ontology.gene_by_name("WDR45").unwrap().to_hpo_set(&ontology),
///     ontology.gene_by_name("TP53").unwrap().to_hpo_set(&ontology),
///     ontology.gene_by_name("CLCN7").unwrap().to_hpo_set(&ontology),
/// ];
///
///
/// let mut cluster = Linkage::union(sets, distance).into_cluster();
/// let first = cluster.next().unwrap();
/// // println!("{:?}", first);
/// // Cluster { idx1: 0, idx2: 3, distance: 0.008127391, size: 2 }
/// assert_eq!(cluster.next().unwrap().len(), 2);
/// assert_eq!(cluster.next().unwrap().len(), 4);
/// assert!(cluster.next().is_none());
/// ```
pub struct Linkage<'a> {
    sets: Vec<Option<HpoSet<'a>>>,
    distance_matrix: DistanceMatrix,
    initial_len: usize,
    clusters: ClusterVec,
}

impl<'a> Linkage<'a> {
    /// Performs union-based hierarchical clustering of `HpoSet`s
    ///
    /// In each iteration, `HpoSet`s are compared to each other based on the
    /// provided `distance` function. `Cluster`s are formed by combining the
    /// 2 closest `HpoSet`s into a single set (forming the union).
    ///
    /// This method becomes exponentially slower with larger lists of sets,
    /// because it merges sets and calculates pairwise similarities
    /// for each term in each set.
    ///
    /// # Examples
    ///
    /// ```rust
    ///
    /// use hpo::Ontology;
    /// use hpo::HpoSet;
    /// use hpo::similarity::GroupSimilarity;
    /// use hpo::utils::Combinations;
    /// use hpo::stats::Linkage;
    ///
    /// // This method can and should utilize parallel processing, e.g.
    /// // using rayon iterators
    /// fn distance(combs: Combinations<HpoSet<'_>>) -> Vec<f32> {
    ///     let sim = GroupSimilarity::default();
    ///     combs.map(|comp| {
    ///         1.0 - sim.calculate(comp.0, comp.1)
    ///     }).collect()
    /// }
    ///
    /// let ontology = Ontology::from_binary("tests/example.hpo").unwrap();
    /// let sets = vec![
    ///     ontology.gene_by_name("KRAS").unwrap().to_hpo_set(&ontology),
    ///     ontology.gene_by_name("WDR45").unwrap().to_hpo_set(&ontology),
    ///     ontology.gene_by_name("TP53").unwrap().to_hpo_set(&ontology),
    ///     ontology.gene_by_name("CLCN7").unwrap().to_hpo_set(&ontology),
    /// ];
    ///
    ///
    /// let mut cluster = Linkage::union(sets, distance).into_cluster();
    /// let first = cluster.next().unwrap();
    /// println!("{:?}", first);
    /// // Cluster { idx1: 0, idx2: 3, distance: 0.16666663, size: 2 }
    /// assert_eq!(cluster.next().unwrap().len(), 2);
    /// assert_eq!(cluster.next().unwrap().len(), 4);
    /// assert!(cluster.next().is_none());
    /// ```
    pub fn union<T, F>(sets: T, distance: F) -> Self
    where
        T: IntoIterator<Item = HpoSet<'a>>,
        F: Fn(Combinations<HpoSet<'_>>) -> Vec<f32>,
    {
        let mut s = Self::new(sets, &distance);
        s.cluster_set_unions(&distance);
        s
    }

    /// Performs single-hierarchical clustering of `HpoSet`s
    ///
    /// `HpoSet`s are compared to each other based on the
    /// provided `distance` function. `Cluster`s are formed by using the minimum
    /// distance of each encompassing set to the comparison set.
    /// This is also known as the Nearest Point Algorithm.
    ///
    /// # Examples
    ///
    /// ```rust
    ///
    /// use hpo::Ontology;
    /// use hpo::HpoSet;
    /// use hpo::similarity::GroupSimilarity;
    /// use hpo::utils::Combinations;
    /// use hpo::stats::Linkage;
    ///
    /// // This method can and should utilize parallel processing, e.g.
    /// // using rayon iterators
    /// fn distance(combs: Combinations<HpoSet<'_>>) -> Vec<f32> {
    ///     let sim = GroupSimilarity::default();
    ///     combs.map(|comp| {
    ///         1.0 - sim.calculate(comp.0, comp.1)
    ///     }).collect()
    /// }
    ///
    /// let ontology = Ontology::from_binary("tests/example.hpo").unwrap();
    /// let sets = vec![
    ///     ontology.gene_by_name("KRAS").unwrap().to_hpo_set(&ontology),
    ///     ontology.gene_by_name("WDR45").unwrap().to_hpo_set(&ontology),
    ///     ontology.gene_by_name("TP53").unwrap().to_hpo_set(&ontology),
    ///     ontology.gene_by_name("CLCN7").unwrap().to_hpo_set(&ontology),
    /// ];
    ///
    ///
    /// let mut cluster = Linkage::single(sets, distance).into_cluster();
    /// let first = cluster.next().unwrap();
    /// println!("{:?}", first);
    /// // Cluster { idx1: 0, idx2: 3, distance: 0.16666663, size: 2 }
    /// assert_eq!(cluster.next().unwrap().len(), 2);
    /// assert_eq!(cluster.next().unwrap().len(), 4);
    /// assert!(cluster.next().is_none());
    /// ```
    pub fn single<T, F>(sets: T, distance: F) -> Self
    where
        T: IntoIterator<Item = HpoSet<'a>>,
        F: Fn(Combinations<HpoSet<'_>>) -> Vec<f32>,
    {
        fn f32_min(v1: Option<&f32>, v2: Option<&f32>) -> f32 {
            if v1.expect("v1 must be `Some`") < v2.expect("v2 must be `Some`") {
                *v1.expect("v1 must be `Some`")
            } else {
                *v2.expect("v2 must be `Some`")
            }
        }

        let mut linkage = Self::new(sets, &distance);
        linkage.arithmetic_cluster(f32_min);
        linkage
    }

    /// Performs complete-hierarchical clustering of `HpoSet`s
    ///
    /// `HpoSet`s are compared to each other based on the
    /// provided `distance` function. `Cluster`s are formed by using the maximum
    /// distance of each encompassing set to the comparison set.
    /// This is also known by the Farthest Point Algorithm or Voor Hees Algorithm.
    pub fn complete<T, F>(sets: T, distance: F) -> Self
    where
        T: IntoIterator<Item = HpoSet<'a>>,
        F: Fn(Combinations<HpoSet<'_>>) -> Vec<f32>,
    {
        fn f32_max(v1: Option<&f32>, v2: Option<&f32>) -> f32 {
            if v1.expect("v1 must be `Some`") > v2.expect("v2 must be `Some`") {
                *v1.expect("v1 must be `Some`")
            } else {
                *v2.expect("v2 must be `Some`")
            }
        }

        let mut linkage = Self::new(sets, &distance);
        linkage.arithmetic_cluster(f32_max);
        linkage
    }

    /// Performs average-hierarchical clustering of `HpoSet`s
    ///
    /// `HpoSet`s are compared to each other based on the
    /// provided `distance` function. `Cluster`s are formed by using the average
    /// distance of both encompassing sets to the comparison set.
    /// This is also called the UPGMA algorithm.
    ///
    /// # Note
    /// This method is not implemented completely correct. Instead of calculating
    /// the average of all distances, it only uses the mean distance of all direct
    /// cluster nodes.
    pub fn average<T, F>(sets: T, distance: F) -> Self
    where
        T: IntoIterator<Item = HpoSet<'a>>,
        F: Fn(Combinations<HpoSet<'_>>) -> Vec<f32>,
    {
        fn mean(v1: Option<&f32>, v2: Option<&f32>) -> f32 {
            (v1.expect("v1 must be `Some`") + v2.expect("v2 must be `Some`")) / 2.0
        }

        let mut linkage = Self::new(sets, &distance);
        linkage.arithmetic_cluster(mean);
        linkage
    }

    /// Returns an Iterator of [`Cluster`] references
    pub fn cluster(&self) -> cluster::Iter {
        self.clusters.iter()
    }

    /// Returns an Iterator of owned [`Cluster`]
    pub fn into_cluster(self) -> cluster::IntoIter {
        self.clusters.into_iter()
    }

    /// Returns the order of the input set items in the final cluster
    ///
    /// # Examples
    ///
    /// ```rust
    ///
    /// use hpo::Ontology;
    /// use hpo::HpoSet;
    /// use hpo::similarity::GroupSimilarity;
    /// use hpo::utils::Combinations;
    /// use hpo::stats::Linkage;
    ///
    /// // This method can and should utilize parallel processing, e.g.
    /// // using rayon iterators
    /// fn distance(combs: Combinations<HpoSet<'_>>) -> Vec<f32> {
    ///     let sim = GroupSimilarity::default();
    ///     combs.map(|comp| {
    ///         1.0 - sim.calculate(comp.0, comp.1)
    ///     }).collect()
    /// }
    ///
    /// let ontology = Ontology::from_binary("tests/example.hpo").unwrap();
    ///
    /// let genes = vec!["KRAS", "WDR45", "ZNRF3", "CLCN7"];
    /// let sets = genes.iter().map(|gene| ontology.gene_by_name(gene).unwrap().to_hpo_set(&ontology));
    ///
    /// let linkage = Linkage::single(sets, distance);
    /// let indicies = linkage.indicies();
    ///
    /// // Similarities:
    /// // CLCN7 WDR45  0.8333334
    /// // CLCN7 ZNRF3  0.66666
    /// // KRAS WDR45   0.5
    /// // KRAS CLCN7   0.4166667
    /// // KRAS ZNRF3   0
    /// // WDR45 ZNRF3  0
    ///
    /// assert_eq!(indicies, vec![1usize, 3usize, 2usize, 0usize]);
    /// for idx in indicies {
    ///    print!("{} ", genes[idx]);
    /// }
    /// // "KRAS ZNRF3 WDR45 CLCN7"
    /// ```
    pub fn indicies(&self) -> Vec<usize> {
        let mut res = Vec::with_capacity(self.initial_len);
        for cluster in &self.clusters {
            if cluster.lhs() < self.initial_len {
                res.push(cluster.lhs());
            }
            if cluster.rhs() < self.initial_len {
                res.push(cluster.rhs());
            }
        }
        res
    }

    fn new<T, F>(sets: T, distance: F) -> Self
    where
        T: IntoIterator<Item = HpoSet<'a>>,
        F: Fn(Combinations<HpoSet<'_>>) -> Vec<f32>,
    {
        let sets: Vec<Option<HpoSet<'a>>> = sets.into_iter().map(Some).collect();
        let len = sets.len();
        let mut s = Self {
            sets,
            distance_matrix: DistanceMatrix::default(),
            initial_len: len,
            clusters: ClusterVec::with_capacity(len),
        };
        // initial calculation of all distances
        s.calculate_initial_distances(distance);
        s
    }

    /// Creates a `DistanceMatrix` with the distances of all `Combinations`
    ///
    /// Note that the `func` return value should be `Distance`, not `Similarity`!
    fn calculate_initial_distances<F: Fn(Combinations<HpoSet<'a>>) -> Vec<f32>>(
        &mut self,
        func: F,
    ) {
        let similarities = func(Combinations::new(&self.sets));

        let index: Vec<Option<usize>> = (0..self.sets.len()).map(Some).collect();
        for ((idx1, idx2), sim) in Combinations::new(&index).zip(similarities.into_iter()) {
            self.distance_matrix.insert((*idx1, *idx2), sim);
        }
    }

    /// Gets the closest two clusters and returns their indicies and their distance
    fn closest_clusters(&self) -> ((usize, usize), f32) {
        self.distance_matrix
            .iter()
            .reduce(|max, elmt| if elmt.1 < max.1 { elmt } else { max })
            .map(|elmt| (*elmt.0, *elmt.1))
            .expect("distance matrix is not empty")
    }

    /// Adds a new cluster
    fn new_cluster(&mut self, key: (usize, usize), dist: f32) {
        self.clusters.push(Cluster::new(
            key.0,
            key.1,
            dist,
            self.size_of_cluster(key.0, key.1),
        ));
    }

    /// Iteratively clusters all sets in the `Linkage` until none are left
    ///
    /// - Finds the 2 clusters/sets with smallest distance
    /// - inactivates them and creates a new, merged set
    /// - calculates the distance between the new set and all other active sets
    /// - appends the `DistanceMatrix` with new distances
    fn cluster_set_unions<F>(&mut self, func: F)
    where
        F: Fn(Combinations<HpoSet<'a>>) -> Vec<f32>,
    {
        loop {
            if self.distance_matrix.is_empty() {
                // All sets are clustered, only one cluster remains
                return;
            }

            // get the indicies of the 2 sets with the smallest distance
            let (key, dist) = self.closest_clusters();

            // create a new cluster with the 2 sets
            self.new_cluster(key, dist);

            // merge the 2 sets and remove them from `self.sets`
            let mut newset = self.sets[key.0]
                .take()
                .expect("set is part of distance matrix and must exist");
            let set2 = self.sets[key.1]
                .take()
                .expect("set is part of distance matrix and must exist");
            newset.extend(&set2);
            self.sets.push(Some(newset));

            // remove all distance scores that include one of the 2 sets
            self.distance_matrix.retain(|(idx1, idx2), _| {
                idx1 != &key.0 && idx1 != &key.1 && idx2 != &key.0 && idx2 != &key.1
            });

            // calculate the distance scores from the new set to all other active sets
            let mut new_combinations = Combinations::new(&self.sets);
            new_combinations.set_to_last();
            let mut distances = func(new_combinations).into_iter();

            // the distance is only calculated for active sets, inactive ones are skipped.
            // due to this, we don't have the proper mapping of index to distance score.
            // to account for this, we're looping through all sets, checking if they are active
            // and then getting the next distance score.
            let last_index = self.sets.len() - 1;
            for (idx, set) in self.sets[..last_index].iter().enumerate() {
                if set.is_some() {
                    self.distance_matrix.insert(
                        (idx, last_index),
                        distances.next().expect("distance score must be present"),
                    );
                }
            }
        }
    }

    /// Iteratively clusters all sets in the `Linkage` until none are left
    ///
    /// - Finds the 2 clusters/sets with smallest distance
    /// - inactivates them and adds a new phantom set at the end
    /// - calculates the distance between the nodes of the new cluster and all
    ///   other clusters based on `func`
    /// - appends the `DistanceMatrix` with new distances
    fn arithmetic_cluster<F>(&mut self, func: F)
    where
        F: Fn(Option<&f32>, Option<&f32>) -> f32,
    {
        loop {
            if self.distance_matrix.is_empty() {
                // All sets are clustered, only one cluster remains
                return;
            }

            // get the indicies of the 2 sets with the smallest distance
            let (key, dist) = self.closest_clusters();

            // create a new cluster with the 2 sets
            self.new_cluster(key, dist);

            // Remove sets
            // For simplicity reasons, add one set to the end as new phantom set
            // (it's pushed to `sets` at the end of the loop)
            let x = self.sets[key.0].take();
            self.sets[key.1].take();

            // add new distance measures to the matrix for
            // each distance between new cluster vs all existing clusters
            // - iterate all existing cluster
            // - compare to new-cluster-1
            // - compare to new-cluster-2
            // - add lower value to distance-matrix
            let new_idx = self.sets.len();
            for (idx, set) in self.sets.iter().enumerate() {
                if idx == key.0 || idx == key.1 {
                    continue;
                }
                if set.is_some() {
                    let distance = match (idx.cmp(&key.0), idx.cmp(&key.1)) {
                        (Ordering::Less, Ordering::Less) => func(
                            self.distance_matrix.get(&(idx, key.0)),
                            self.distance_matrix.get(&(idx, key.1)),
                        ),
                        (Ordering::Less, Ordering::Greater) => func(
                            self.distance_matrix.get(&(idx, key.0)),
                            self.distance_matrix.get(&(key.1, idx)),
                        ),
                        (Ordering::Greater, Ordering::Less) => func(
                            self.distance_matrix.get(&(key.0, idx)),
                            self.distance_matrix.get(&(idx, key.1)),
                        ),
                        (Ordering::Greater, Ordering::Greater) => func(
                            self.distance_matrix.get(&(key.0, idx)),
                            self.distance_matrix.get(&(key.1, idx)),
                        ),
                        (Ordering::Equal, _) | (_, Ordering::Equal) => {
                            unreachable!("Cannot be reached")
                        }
                    };
                    self.distance_matrix.insert((idx, new_idx), distance);
                }
            }
            // remove all distance scores that include one of the 2 sets
            self.distance_matrix.retain(|(idx1, idx2), _| {
                idx1 != &key.0 && idx1 != &key.1 && idx2 != &key.0 && idx2 != &key.1
            });

            self.sets.push(x);
        }
    }

    /// Returns the size of the cluster, the sum of sizes of both nodes
    fn size_of_cluster(&self, idx1: usize, idx2: usize) -> usize {
        (if idx1 < self.initial_len {
            1
        } else {
            self.clusters
                .get(idx1 - self.initial_len)
                .expect("idx is guaranteed to be in cluster")
                .len()
        }) + (if idx2 < self.initial_len {
            1
        } else {
            self.clusters
                .get(idx2 - self.initial_len)
                .expect("idx is guaranteed to be in cluster")
                .len()
        })
    }

    /// Returns an Iterator of [`Cluster`] references
    pub fn iter(&self) -> cluster::Iter {
        self.cluster()
    }
}

impl<'a> IntoIterator for &'a Linkage<'a> {
    type Item = &'a Cluster;
    type IntoIter = cluster::Iter<'a>;
    fn into_iter(self) -> Self::IntoIter {
        self.iter()
    }
}

impl IntoIterator for Linkage<'_> {
    type Item = Cluster;
    type IntoIter = cluster::IntoIter;
    fn into_iter(self) -> Self::IntoIter {
        self.into_cluster()
    }
}
