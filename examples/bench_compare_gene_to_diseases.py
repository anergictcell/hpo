from pyhpo import Ontology, HPOSet
from pyhpo.annotations import Gene

_ = Ontology()

count = 0

gba = Gene.get("GBA1")
gba_set = HPOSet.from_queries(gba.hpo)

for disease in Ontology.omim_diseases:
    s = HPOSet.from_queries(disease.hpo)
    if gba_set.similarity(
        s,
        kind="omim",
        method="graphic",
        combine="funSimAvg"
    ) > 0.6:
        count += 1

print(f"Highly similar: {count}")
