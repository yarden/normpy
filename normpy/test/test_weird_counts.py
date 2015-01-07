##
## Test a weird counts dataset
##
import os
import pandas

import normpy
import normpy.utils as utils
import normpy.experiment as experiment
import normpy.normalizers as normalizers

import matplotlib.pylab as plt

def test_weird_counts():
    samples = {"counts_A": "counts_A",
               "counts_B": "counts_B",
               "counts_C": "counts_C",
               "counts_D": "counts_D",
               "counts_E": "counts_E"}
    counts = pandas.read_table("./data/weird_counts.txt", sep="\t")
    counts = counts.set_index("gene_id")
    exp_obj = experiment.Experiment(counts, samples)
    # DESeq normalization
    deseq_norm = normalizers.norm_deseq(exp_obj)
    # TMM normalization
    tmm_norm = normalizers.norm_tmm(exp_obj)
    # Compare results for one gene
    g = "ENSMUSG00000020140"
    print "DESeq results for %s" %(g)
    print deseq_norm.ix[g]
    print "=" * 10
    print "TMM resutls for %s" %(g)
    print tmm_norm.ix[g]
    # fold difference in counts for sample C
    print "Fold difference in sample C (TMM / DESeq): %.4f" \
          %(tmm_norm.ix[g]["counts_C"] / deseq_norm.ix[g]["counts_C"])
    # fold relative difference in counts, comparing A to C
    tmm_A_vs_C = tmm_norm.ix[g]["counts_A"] / tmm_norm.ix[g]["counts_C"]
    deseq_A_vs_C = deseq_norm.ix[g]["counts_A"] / deseq_norm.ix[g]["counts_C"]
    print "TMM fold change A vs. C: %.2f" %(tmm_A_vs_C)
    print "DESeq fold change A vs. C: %.2f" %(deseq_A_vs_C)
    

def main():
    test_weird_counts()



if __name__ == "__main__":
    main()
