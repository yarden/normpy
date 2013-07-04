##
## Plotting utilities
##
import os
import scipy
import numpy as np

import pandas

import normpy
import normpy.utils as utils

import matplotlib
import matplotlib.pylab as plt


def plot_fcs(normed_df, unnormed_df, pair, basename):
    """
    Plot fold changes for normed and unnormed dataframes.

    Parameters:
    -----------
    normed_df: Normalized dataframe of values
    unnormed_df: Unnormalized dataframe of values
    pair: Tuple containing the two columns to use
    to compute fold change. Fold change is first
    sample divided by second.
    """
    if (pair[0] not in normed_df.columns) or \
       (pair[1] not in normed_df.columns):
        raise Exception, "One of the pairs is not in normed df."
    if (pair[0] not in unnormed_df.columns) or \
       (pair[1] not in unnormed_df.columns):
        raise Exception, "One of the pairs is not in unnormed.df"
    normed_fc = \
        np.log2(normed_df[pair[0]]) - np.log2(normed_df[pair[1]])
    unnormed_fc = \
        np.log2(unnormed_df[pair[0]]) - np.log2(unnormed_df[pair[1]])
    fc_df = pandas.DataFrame({"normed_fc": normed_fc,
                              "unnormed_fc": unnormed_fc})
    # Remove nans/infs etc.
    pandas.set_option('use_inf_as_null', True)
    fc_df = fc_df.dropna(how="any", subset=["normed_fc",
                                            "unnormed_fc"])
    plt.figure()
    fc_df.hist(color="k", bins=40)
    plt.suptitle("%s vs. %s" %(pair[0], pair[1]))
    plt.xlabel("Fold change (log2)")
    save_fig(basename)
    pandas.set_option('use_inf_as_null', False)    


def save_fig(basename, ext="pdf"):
    """
    Save figure to plots dir.
    """
    utils.make_dir(utils.get_plots_dir())
    plot_fname = os.path.join(utils.get_plots_dir(),
                              "%s.%s" %(basename, ext))
    plt.savefig(plot_fname)
