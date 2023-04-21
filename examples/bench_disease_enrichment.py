from pyhpo import Ontology
from pyhpo.annotations import Gene
from pyhpo.set import HPOSet
from pyhpo.stats import EnrichmentModel


GENES = ["GBA1", "NPC1", "EZH2", "DMD", "MUC7", "ARID1B"]

_ = Ontology()
enrichment = EnrichmentModel('omim')

count = 0
for gene in GENES:
    ci = HPOSet.from_queries(Gene.get(gene).hpo)
    for res in enrichment.enrichment(method='hypergeom', hposet=ci):
        if res["enrichment"] < 0.0000005:
            count += 1

print(f"Highly enriched: {count}")
