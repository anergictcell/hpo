# Release checklist

This document contains the workflows to follow for all changes and releases to `hpo`. 
The worklow assures that the `main` branch always holds a functional version of `hpo` with all tests passing. The `main` branch can be ahead of the official `crates.io` release. New versions for `crates.io` releases are created independently of the regular updates and will contain all changes present in the `main` branch at that point. My goal is to automate the version bump and release process using Github Actions at some point.

This procedure is just a suggestion at this point and can be modified if needs arise.


## Regular updates / Normal development

- [ ] Develop in a dedicated branch (or your own fork): `git checkout -b <MY_FEATURE_NAME>`
- [ ] Rebase onto `main`: `git rebase main <MY_FEATURE_NAME>`
- [ ] Double check for good code, sensible API and well-explained docs
- [ ] Run format, clippy, tests and doc-generation: `cargo fmt --check && cargo clippy && cargo test && cargo doc`
- [ ] Push to remote: `git push -u origin <MY_FEATURE_NAME>`
- [ ] Create merge/pull request to `main` branch
- [ ] Once CICD passes, changes are merged to `main`


## Version bumps

- [ ] Make dedicated branch named after version: `git checkout main && git pull && git checkout -b release/<MAJOR>.<MINOR>.<PATCH>`
- [ ] Update Cargo.toml with new version
- [ ] Update dependencies if needed and possible
- [ ] Check if README or docs need update
- [ ] Add Changelog summary of changes
- [ ] Run format, clippy, tests and doc-generation: `cargo fmt --check && cargo clippy && cargo test && cargo doc`
- [ ] add git tag with version: `git tag v<MAJOR>.<MINOR>.<PATCH>`
- [ ] push to remote, also push tags: `git push -u origin release/<MAJOR>.<MINOR>.<PATCH> && git push --tags`
- [ ] Merge into main
- [ ] update main branch locally: `git checkout main && git pull`
- [ ] release to cargo: `cargo publish`
