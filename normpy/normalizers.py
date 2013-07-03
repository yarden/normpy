##
## Normalizers
##
import numpy as np
import scipy

# Wrapper functions for normalizers
def norm_uq():
    """
    Upper quantile.
    """
    pass


def norm_tmm():
    """
    TMM.
    """
    pass


def norm_tc(exp_obj):
    """
    Normalize by total counts. Return new counts dataframe.

    Parameters:
    -----------
    - exp_obj: experiment object. Normalizes by sum of total counts
    or by library size if given.
    """
    # Copy the dataframe (maybe a better way to do this...?)
    norm_counts_df = exp_obj.counts_df.copy()
    # Mapping from column to the normalizing factor
    total_counts = {}
    if exp_obj.lib_sizes is not None:
        # Use the library sizes to normalize
        for sample_col in exp_obj.lib_sizes:
            total_counts[sample_col] = exp_obj.lib_sizes[sample_col]
    else:
        # Otherwise, use the sum of the columns
        for sample_col in exp_obj.counts_df.columns:
            total_counts[sample_col] = exp_obj.counts_df[sample_col].sum()
    # Normalize
    for sample_col in norm_counts_df.columns:
        norm_denom = float(total_counts[sample_col])
        if norm_denom == 0:
            raise Exception, "Sample %s is an empty library!" %(sample_col)
        norm_counts_df[sample_col] = norm_counts_df[sample_col] / norm_denom
    return norm_counts_df
        

def norm_lowess():
    """
    Lowess.
    """
    pass


def norm_rpkm():
    """
    RPKM [optional].
    """
    pass


# Methods:
#  UQ (Upper Quartile)
#  TMM (Trimmed Mean of M values)
#  TC (Total counts)
#  LOWESS
#  [optional] RPKM (RPKM)

