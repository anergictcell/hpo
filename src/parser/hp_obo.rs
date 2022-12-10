// //! This will work the following way:
// //! 1. Add new terms
// //! 2. Replace obsolete ones
// //! 3. Connect parents


/*
Different options provided by lib:
1. Read from OBO
2. Have 
    - one file with HPOTerm(id,name)
    - one file with term->parents mapping

*/

// // }
use log::{trace, warn};

use crate::parser::Path;
use std::fs;

use crate::{HpoTermId, term::internal::HpoTermInternal, Ontology};

type Connections = Vec<(HpoTermId, HpoTermId)>;

pub (super) fn read_obo_file<P: AsRef<Path>>(filename: P, ontology: &mut Ontology)  {
    // stores tuples of Term - Parent
    let mut connections: Connections = Vec::new();

    let file_content = fs::read_to_string(filename).unwrap();
    for term in file_content.split("\n\n") {
        if let Some(term) = term.strip_prefix("[Term]\n") {
            if let Some(raw_term) = term_from_obo(term){
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
}

fn term_from_obo(term: &str) -> Option<HpoTermInternal> {
    let mut id: Option<&str> = None;
    let mut name: Option<&str> = None;
    for line in term.lines() {
        match parse_line(line) {
            ("id", value) => id = Some(value),
            ("name", value) => name = Some(value),
            _ => ()
        }
        if let (Some(id), Some(name)) = (id, name) {
            return Some(HpoTermInternal::try_new(id, name).unwrap())
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
        let mut ont = Ontology::empty();
        read_obo_file("tests/small.obo", &mut ont);

        assert_eq!(ont.len(), 4);

        assert!(ont.hpo(&1u32.into()).is_none());

        assert_eq!(ont.hpo(&219u32.into()).unwrap().parents().count(), 2);

        assert_eq!(ont.hpo(&218u32.into()).unwrap().parents().count(), 1);

        assert_eq!(ont.hpo(&217u32.into()).unwrap().parents().count(), 0);

        assert_eq!(ont.hpo(&217u32.into()).unwrap().children().count(), 2);
    }
}