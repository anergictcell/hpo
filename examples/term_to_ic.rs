use hpo::Ontology;

fn main() {
    let ontology = Ontology::from_standard("./example_data/").unwrap();

    for term in ontology.hpos() {
        println!(
            "{}\t{}\t{}",
            term.id(),
            term.information_content().gene(),
            term.information_content().omim_disease(),
        );
    }
}

/*
Python comparison code
from pyhpo import Ontology

_ = Ontology()

with open("term2ic.py.txt", "w") as fh:
    for term in Ontology:
        _ = fh.write(f"{term.id}\t{term.information_content.gene}\t{term.information_content.omim}\n")

*/

/*
Diff'ing:

diff \
<(awk '{print $1, substr($2,1,5), substr($3,1,5)}' example_data/term2ic.rs.txt) \
<(awk '{print $1, substr($2,1,5), substr($3,1,5)}' example_data/term2ic.py.txt)

*/