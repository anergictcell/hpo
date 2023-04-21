use criterion::{black_box, criterion_group, criterion_main, Criterion};

use hpo::HpoTermId;
use hpo::Ontology;

fn ancestors(ontology: &Ontology, times: usize) -> usize {
    let mut count = 0;
    let mut terms: (HpoTermId, HpoTermId) = (
        HpoTermId::try_from("HP:0000001").unwrap(),
        HpoTermId::try_from("HP:0000001").unwrap(),
    );
    ontology
        .into_iter()
        .skip(800)
        .take(times)
        .for_each(|term1| {
            for term2 in ontology.into_iter().skip(500).take(times) {
                let overlap = term1.common_ancestor_ids(&term2).len();
                if overlap > count {
                    count = overlap;
                    terms = (term1.id(), term2.id());
                }
            }
        });
    count
}

fn union_ancestors(ontology: &Ontology, times: usize) -> usize {
    let mut count = 0;
    let mut terms: (HpoTermId, HpoTermId) = (
        HpoTermId::try_from("HP:0000001").unwrap(),
        HpoTermId::try_from("HP:0000001").unwrap(),
    );
    ontology
        .into_iter()
        .skip(800)
        .take(times)
        .for_each(|term1| {
            for term2 in ontology.into_iter().skip(500).take(times) {
                let overlap = term1.union_ancestor_ids(&term2).len();
                if overlap > count {
                    count = overlap;
                    terms = (term1.id(), term2.id());
                }
            }
        });
    count
}

fn ancestors_benchmark(c: &mut Criterion) {
    let ontology = Ontology::from_binary("tests/ontology.hpo").unwrap();

    c.bench_function("common-ancestors 500", |b| {
        b.iter(|| ancestors(black_box(&ontology), black_box(500)))
    });
}

fn union_ancestors_benchmark(c: &mut Criterion) {
    let ontology = Ontology::from_binary("tests/ontology.hpo").unwrap();

    c.bench_function("union-ancestors 500", |b| {
        b.iter(|| union_ancestors(black_box(&ontology), black_box(500)))
    });
}

criterion_group!(
    common_ancestors,
    ancestors_benchmark,
    union_ancestors_benchmark
);
criterion_main!(common_ancestors);
