use criterion::{black_box, criterion_group, criterion_main, Criterion};
use rayon::prelude::*;

use hpo::term::InformationContentKind;

use hpo::similarity::Builtins;
use hpo::Ontology;

fn bench(ontology: &Ontology, times: usize) -> usize {
    let ic = Builtins::GraphIc(InformationContentKind::Omim);
    let mut count = 0usize;
    for term1 in ontology.into_iter().skip(1000).take(times) {
        for term2 in ontology.into_iter().skip(500).take(times) {
            let overlap = term1.similarity_score(&term2, &ic);
            if overlap > 0.7 {
                count += 1;
            }
        }
    }
    count
}

fn parallel(ontology: &Ontology, times: usize) -> usize {
    let ic = Builtins::GraphIc(InformationContentKind::Omim);

    ontology
        .into_iter()
        .skip(1000)
        .take(times)
        .par_bridge()
        .map(|term1| {
            let mut count = 0usize;
            for term2 in ontology.into_iter().skip(1000).take(times) {
                let overlap = term1.similarity_score(&term2, &ic);
                if overlap > 0.7 {
                    count += 1;
                }
            }
            count
        })
        .sum()
}

fn graphic_benchmark(c: &mut Criterion) {
    let ontology = Ontology::from_standard("./example_data/").unwrap();

    c.bench_function("graphic 100", |b| {
        b.iter(|| bench(black_box(&ontology), black_box(100)))
    });

    c.bench_function("graphic-parallel 1000", |b| {
        b.iter(|| parallel(black_box(&ontology), black_box(1000)))
    });
}

criterion_group!(similarity, graphic_benchmark);
criterion_main!(similarity);
