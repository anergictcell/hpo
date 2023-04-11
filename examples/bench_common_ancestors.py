from pyhpo import Ontology

_ = Ontology()

common = 0
terms = (None, None)
for term1 in Ontology:
    for term2 in list(Ontology)[0:10000]:
        overlap = term1.common_ancestors(term2)
        if len(overlap) > common:
            common = len(overlap)
            terms = (term1, term2)

print("Terms {} and {} have {} overlaps".format(
    terms[0].id,
    terms[1].id,
    common
))
