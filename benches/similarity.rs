use criterion::{black_box, criterion_group, criterion_main, Criterion};
use rayon::prelude::*;

use hpo::term::InformationContentKind;

use hpo::similarity::Builtins;
use hpo::Ontology;

fn graphic_sequential(ontology: &Ontology, times: usize) -> usize {
    let ic = Builtins::GraphIc(InformationContentKind::Omim);
    let mut count = 0usize;
    for term1 in ontology.into_iter().skip(300).take(times) {
        for term2 in ontology.into_iter().skip(500).take(times) {
            let overlap = term1.similarity_score(&term2, &ic);
            if overlap > 0.7 {
                count += 1;
            }
        }
    }
    count
}

fn graphic_parallel(ontology: &Ontology, times: usize) -> usize {
    let ic = Builtins::GraphIc(InformationContentKind::Omim);

    ontology
        .into_iter()
        .skip(800)
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
    let ontology = Ontology::from_binary("tests/ontology.hpo").unwrap();

    c.bench_function("graphic 300", |b| {
        b.iter(|| graphic_sequential(black_box(&ontology), black_box(300)))
    });

    c.bench_function("graphic-parallel 1000", |b| {
        b.iter(|| graphic_parallel(black_box(&ontology), black_box(1000)))
    });
}


criterion_group!(similarity, graphic_benchmark);
criterion_main!(similarity);
