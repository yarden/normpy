##
## Unit tests
##
import os

import pandas

import numpy as np

import normpy
import normpy.utils as utils
import normpy.edgeR_utils as edgeR_utils
import normpy.experiment as experiment
import normpy.normalizers as normalizers
import normpy.plot_utils as plot_utils

import matplotlib
import matplotlib.pylab as plt

import scipy
import scipy.stats

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


def test_deseq():
    """
    Calls DESeq normalization. Prints raw counts then normed counts.
    """
    counts_fname = utils.load_testdata("pasilla")
    # Consider only a subset of the samples
    samples = OrderedDict()
    samples["Untreated 1"] = "untreated1"
    samples["Untreated 2"] = "untreated2"
    exp_obj = experiment.Experiment(counts_fname, samples)
    norm_counts_df = normalizers.norm_deseq(exp_obj)
    print "\nDESeq Testing:"
    print "--------------"
    print "Pre-normalized counts: "
    print exp_obj.counts_df.head()
    print "Normalized counts: "
    print norm_counts_df.head()


def test_lowess():
    """
    Tests lowess normalization. 
    """
    counts_fname = utils.load_testdata("pasilla")
    # Consider only a subset of the samples
    samples = OrderedDict()
    samples["Untreated 1"] = "untreated1"
    samples["Untreated 2"] = "untreated2"
    exp_obj = experiment.Experiment(counts_fname, samples)
    pairs = [["untreated1", "untreated2"]]
    norm_df, unnorm_df = normalizers.norm_ma_lowess(exp_obj, pairs)
    print "\nLowess Testing:"
    print "--------------"
    print "Pre-normalized values: "
    print unnorm_df.head()
    print "Normalized counts: "
    print norm_df.head()
    # Compare LOWESS normalized to total counts
    pair = ["untreated1", "untreated2"]
    plot_utils.plot_fcs(norm_df, unnorm_df, pair, "lowess_test")


def R_test_lowess():
    """
    Tests lowess normalization. 
    """
    counts_fname = utils.load_testdata("pasilla")
    # Consider only a subset of the samples
    samples = OrderedDict()
    samples["Untreated 1"] = "untreated1"
    samples["Untreated 2"] = "untreated2"
    exp_obj = experiment.Experiment(counts_fname, samples)
    pairs = [["untreated1", "untreated2"]]
    norm_df, unnorm_df = normalizers.R_norm_ma_lowess(exp_obj, pairs)
    print "\nLowess Testing:"
    print "--------------"
    print "Pre-normalized values: "
    print unnorm_df.head()
    print "Normalized counts: "
    print norm_df.head()
    # Compare LOWESS normalized to total counts
    pair = ["untreated1", "untreated2"]
    plot_utils.plot_fcs(norm_df, unnorm_df, pair, "R_lowess_test")


def test_tmm():
    """
    Calls TMM normalization. Prints raw counts then normed counts.
    """
    counts_fname = utils.load_testdata("pasilla")
    # Consider only a subset of the samples
    samples = OrderedDict()
    samples["Untreated 1"] = "untreated1"
    samples["Untreated 2"] = "untreated2"
    exp_obj = experiment.Experiment(counts_fname, samples)
    norm_counts_df = normalizers.norm_tmm(exp_obj)
    print "\nTMM Testing:"
    print "--------------"
    print "Pre-normalized counts: "
    print exp_obj.counts_df.head()
    print "Normalized counts: "
    print norm_counts_df.head()


def test_quantile():
    """
    Calls quantile normalization. Prints raw counts then normed counts.
    """
    counts_fname = utils.load_testdata("pasilla")
    # Consider only a subset of the samples
    samples = OrderedDict()
    samples["Untreated 1"] = "untreated1"
    samples["Untreated 2"] = "untreated2"
    exp_obj = experiment.Experiment(counts_fname, samples)
    norm_counts_df = normalizers.norm_q(exp_obj)
    print "\nQuantile Testing:"
    print "--------------"
    print "Pre-normalized counts: "
    print exp_obj.counts_df.head()
    print "Normalized counts: "
    print norm_counts_df.head()


def test_quantile_vs_tmm():
    """
    Test quantile normalization versus TMM
    in rank correlation of genes.
    """
    counts_fname = utils.load_testdata("pasilla")
    # Consider only a subset of the samples
    samples = OrderedDict()
    samples["Untreated 1"] = "untreated1"
    samples["Untreated 2"] = "untreated2"
    exp_obj = experiment.Experiment(counts_fname, samples)
    quantile_counts_df = normalizers.norm_q(exp_obj)
    tmm_counts_df = normalizers.norm_tmm(exp_obj)
    print "\nQuantile versus TMM Testing:"
    print "--------------"
    print "Normalized quantile counts: "
    print quantile_counts_df.head()
    print "Normalized TMM counts: "
    print tmm_counts_df.head()
    print "Correlating the genes."
    # Merge the dataframes together, indexing by gene
    combined_df = pandas.merge(quantile_counts_df, tmm_counts_df,
                               left_index=True,
                               right_index=True,
                               suffixes=["_q", "_tmm"],
                               how="outer")
    # Get log of counts: get rid of infinite values
    log_counts_df = combined_df.apply(np.log2).replace([-np.inf, np.inf],
                                                       np.nan)
    print "Combined dataframe: "
    print combined_df.head()
    print "Combined log dataframe: "
    print log_counts_df.head()
    # Plot correlation
    from pandas.tools.plotting import scatter_matrix
    scatter_matrix(log_counts_df, alpha=0.2, figsize=(8, 7))
    plot_utils.save_fig("quantile_vs_tmm_corr", ext="png")

    
def main():
    test_edgeR()
    test_tc()
    test_deseq()
    #test_lowess()
    R_test_lowess()
    test_tmm()
    test_quantile()
    test_quantile_vs_tmm()


if __name__ == "__main__":
    main()

    
