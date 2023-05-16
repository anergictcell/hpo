import itertools
import sys

from pyhpo import Ontology

_ = Ontology()

count = 0

combinations = itertools.product(list(Ontology), list(Ontology)[0:1000])

if len(sys.argv) == 1:
    for c in combinations:
        if c[0].similarity_score(c[1], kind="omim", method="graphic") > 0.9:
            count += 1
else:
    try:
        from pyhpo.helper import batch_similarity
    except ModuleNotFoundError as err:
        print("This is only available when using hpo3")
        raise err
    count = len([x for x in batch_similarity(list(combinations), kind="omim", method="graphic") if x > 0.9])


print(f"Comparisons above 0.9 {count}")
