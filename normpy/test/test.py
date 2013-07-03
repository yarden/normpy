##
## Unit tests
##
import os
import normpy
import normpy.utils as utils
import normpy.edgeR_utils as edgeR_utils
import normpy.experiment as experiment
import normpy.normalizers as normalizers

from collections import OrderedDict

def test_edgeR():
    """
    Test call to EdgeR as a reference.
    """
    #counts_fname = "/home/yarden/jaen/normpy/normpy/test/data/pasilla_gene_counts.tsv"
    counts_fname = utils.load_testdata("pasilla")
    # Consider only a subset of the samples
    samples = OrderedDict()
    samples["Untreated 1"] = "untreated1"
    samples["Untreated 2"] = "untreated2"
    # Make an experiment object out of them
    exp_obj = experiment.Experiment(counts_fname, samples)
    edgeR_obj = edgeR_utils.EdgeR(exp_obj)
    edgeR_obj.norm_expr_vals(ref_col="untreated1")


def test_tc():
    """
    Test total counts normalization.
    """
    counts_fname = utils.load_testdata("pasilla")
    # Consider only a subset of the samples
    samples = OrderedDict()
    samples["Untreated 1"] = "untreated1"
    samples["Untreated 2"] = "untreated2"
    # Make an experiment object out of them
    exp_obj = experiment.Experiment(counts_fname, samples)
    # Normalization without explicit library sizes
    norm_counts_df = normalizers.norm_tc(exp_obj)
    print "Pre-normalized counts: "
    print exp_obj.counts_df.head()
    print "Normalized counts: "
    print norm_counts_df.head()
    # Normalization with library sizes
    exp_obj.lib_sizes = {"untreated1": 100,
                         "untreated2": 200}
    norm_counts_df = normalizers.norm_tc(exp_obj)
    print "Normalized with library size: "
    print norm_counts_df.head()


def main():
    test_edgeR()
    test_tc()


if __name__ == "__main__":
    main()

    
