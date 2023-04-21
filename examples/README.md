# HPO - Examples
This folder contains some examples to demonstrate how to use the `hpo` library.
Some of those examples were used during the development to verify the correctness of the business logic. Some examples write huge amounts of output to the console, so check before you run them.

The examples are using different input source data and CLI arguments. Most rely on the source data format provided by JAX, saved to a `example_data` subfolders.

Others require a `ontology.hpo` binary data file. The binary data can be generated manually using the `obo_to_bin` example:

```bash
cargo run --release --example obo_to_bin <PATH TO JAX SOURCE DATA> <OUTPUT_FILE>

# e.g.:
cargo run --release --example obo_to_bin example_data/ example_data/ontology.hpo
```
(There is also an `example.hpo` binary file available in the `tests` subfolder)

Some examples contain the corresponding Python code for `PyHPO` in the comments to allow easier reproducability and show expected runtimes.


# Benchmarks
All examples that are named `bench_*` are used for benchmarking against PyHPO. Each has a corresponding Python script as well.

## Run the Rust benchmarks
```bash
# Build each example binary once before the benchmark
for file in examples/bench_*.rs;
do
    bench=$(basename $file .rs)
    cargo run --release --example ${bench} example_data/ontology.hpo parallel
done

# Benchmark the prebuild binaries
for file in examples/bench_*.rs;
do
    bench=$(basename $file .rs)
    echo ${bench}
    /usr/bin/time target/release/examples/${bench} example_data/ontology.hpo && \
    /usr/bin/time target/release/examples/${bench} example_data/ontology.hpo parallel
done
```

## Run the Python benchmarks
```bash
for file in examples/bench_*.py;
do
    echo ${file}
    /usr/bin/time python ${file}
done
```