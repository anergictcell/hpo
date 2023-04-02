use hpo::Ontology;
use std::io::Read;
use std::{fs::File, time::Duration};

use criterion::{black_box, criterion_group, criterion_main, Criterion};

fn build_ontology_benchmark(c: &mut Criterion) {
    let mut bytes = Vec::new();
    File::open("tests/ontology.hpo")
        .expect("ontology.hpo cannot be opened")
        .read_to_end(&mut bytes)
        .expect("cannot read ontology.hpo");
    c.bench_function("build ontology", |b| {
        b.iter(|| {
            Ontology::from_bytes(black_box(&bytes[..]))
                .expect("requires valid bytes")
                .len()
        })
    });
}

criterion_group! {
    name = ontology;
    config = Criterion::default().sample_size(20).measurement_time(Duration::from_secs(10));
    targets = build_ontology_benchmark
}
criterion_main!(ontology);
