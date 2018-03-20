"""
The purpose of this script is to randomly sample within each subtype to get
an id for the min identity within a given subtype.

SKH 20180320
"""

import pandas as pd
from Bio import SeqIO
from collections import Counter
import itertools
import glob
from pymodules.utils import *
import numpy as np
import scipy
# for pair in itertools.combinations(group["seq_id"], 2)

# Set up for translation and pairwise
# codon to index, add codon `---`
max_codon_index = (max(list(constants.CODON_TO_INDEX.values())))
if "---" not in list(constants.CODON_TO_INDEX.keys()):
    constants.CODON_TO_INDEX["---"] = max_codon_index + 1
# index to AA, add AA `-`
max_AA_index = max(list(constants.INDEX_TO_AA.keys()))
if "-" not in list(constants.INDEX_TO_AA.values()):
    constants.INDEX_TO_AA[max_AA_index+1] = "-"
# condon to AA, add mapping of `---` to `-`
constants.CODON_TO_AA = scipy.append(constants.CODON_TO_AA, [[max_AA_index+1]])

def main():
    np.random.seed(0)
    n_samples = 50000
    final = {"HA":[], "min_aa_id":[]}
    NUMBERS = [x for x in range(1, 19) if x not in [15, 17, 18]]
    for HA in NUMBERS:
        m = 1.0
        seqs = SeqIO.to_dict(SeqIO.parse("aligned/HA_{0}_hybrid_align.fasta".format(HA), "fasta"))
        ids = list(seqs.keys())
        ids.sort()
        for i in range(n_samples):
            sample = np.random.choice(ids, 2, replace=False)
            p1 = translate_with_gaps(str(seqs[sample[0]].seq))
            p2 = translate_with_gaps(str(seqs[sample[1]].seq))
            aa_id = [1 if p1[x] == p2[x] else 0 for
                          x in range(len(p1))]
            aa_id = sum(aa_id)/len(aa_id)
            if aa_id < m:
                m = aa_id
        final["HA"].append(HA)
        final["min_aa_id"].append(m)
        print(HA, m)
    final = pd.DataFrame(final)
    final.to_csv("_temp_subtype_aa_id.csv", index=False)




if __name__ == '__main__':
    main()
