from pyhpo import Ontology

_ = Ontology()

count = 0

for term_1 in Ontology:
    for term_2 in list(Ontology)[0:1000]:
        if term_1.similarity_score(term_2, kind="omim", method="graphic") > 0.9:
            count += 1

print(f"Comparisons above 0.9 {count}")
