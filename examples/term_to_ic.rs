//! Prints every term and its associated InformationContent

use hpo::Ontology;

fn main() {
    let ontology = Ontology::from_binary("tests/ontology.hpo").unwrap();

    for term in ontology.hpos() {
        println!(
            "{}\t{}\t{}\t{}",
            term.id(),
            term.information_content().gene(),
            term.information_content().omim_disease(),
            term.information_content().orpha_disease(),
        );
    }
}

/*
Python comparison code
from pyhpo import Ontology

_ = Ontology()

with open("term2ic.py.txt", "w") as fh:
    for term in Ontology:
        _ = fh.write(f"{term.id}\t{term.information_content.gene}\t{term.information_content.omim}\t{term.information_content.orpha}\n")

*/

/*
Diff'ing (only use the first 2 decimal digits to ignore rounding errors:

diff \
<(sort example_data/term2ic.rs.txt | awk '{print $1, substr($2,1,4), substr($3,1,4), substr($4,1,4)}') \
<(sort example_data/term2ic.py.txt | awk '{print $1, substr($2,1,4), substr($3,1,4), substr($4,1,4)}')

*/
