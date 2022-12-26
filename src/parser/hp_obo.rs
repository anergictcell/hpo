use log::{trace, warn};

use crate::{parser::Path, HpoError, HpoResult};
use std::fs;

use crate::{term::internal::HpoTermInternal, HpoTermId, Ontology};

/// Links terms to each other (Child - Parent)
///
/// This type is used to store the Term - Parent dependencies
/// for all HPO-Terms
type Connections = Vec<(HpoTermId, HpoTermId)>;

/// Parses the `hp.obo` file as provided by Jax
///
/// It extracts the `HpoTermId`, the name and the list of parents
/// and inserts them into the [`Ontology`]
///
/// Once the parsing is finshed, the indirect parents for each term
/// are cached.
///
/// If you use this function you cannot add additional terms or
/// parents afterwards, since all dependency data will be already cached.
pub(super) fn read_obo_file<P: AsRef<Path>>(filename: P, ontology: &mut Ontology) -> HpoResult<()> {
    // stores tuples of Term - Parent
    let mut connections: Connections = Vec::new();

    let file_content = match fs::read_to_string(&filename) {
        Ok(content) => content,
        Err(_) => {
            return Err(HpoError::CannotOpenFile(
                filename.as_ref().display().to_string(),
            ))
        }
    };
    for term in file_content.split("\n\n") {
        if let Some(term) = term.strip_prefix("[Term]\n") {
            if let Some(raw_term) = term_from_obo(term) {
                let id = ontology.add_term(raw_term);
                add_connections(&mut connections, term, id);
            } else {
                warn!("Unable to parse: {}", term);
            }
        } else {
            trace!("Ignoring: {}", term);
        }
    }

    for (child, parent) in connections {
        ontology.add_parent(parent, child);
    }

    ontology.create_cache();
    Ok(())
}

fn term_from_obo(term: &str) -> Option<HpoTermInternal> {
    let mut id: Option<&str> = None;
    let mut name: Option<&str> = None;
    for line in term.lines() {
        match parse_line(line) {
            ("id", value) => id = Some(value),
            ("name", value) => name = Some(value),
            _ => (),
        }
        if let (Some(id), Some(name)) = (id, name) {
            return Some(HpoTermInternal::try_new(id, name).unwrap());
        }
    }
    None
}

fn add_connections(connections: &mut Connections, term: &str, id: HpoTermId) {
    for line in term.lines() {
        if let Some(value) = line.strip_prefix("is_a: ") {
            if let Some((term_id, _)) = value.split_once(' ') {
                connections.push((id, HpoTermId::try_from(term_id).unwrap()))
            } else {
                println!("Unable to parse HPO ID from {}", value)
            }
        }
    }
}

fn parse_line(line: &str) -> (&str, &str) {
    line.split_once(": ").expect("unable to parse line")
}

#[cfg(test)]
mod test {
    // use std::fs;

    use super::*;
    use crate::Ontology;

    #[test]
    fn split_terms() {
        let mut ont = Ontology::default();
        read_obo_file("tests/small.obo", &mut ont).unwrap();

        assert_eq!(ont.len(), 4);

        assert!(ont.hpo(&1u32.into()).is_none());

        assert_eq!(ont.hpo(&219u32.into()).unwrap().parents().count(), 2);

        assert_eq!(ont.hpo(&218u32.into()).unwrap().parents().count(), 1);

        assert_eq!(ont.hpo(&217u32.into()).unwrap().parents().count(), 0);

        assert_eq!(ont.hpo(&217u32.into()).unwrap().children().count(), 2);
    }
}
