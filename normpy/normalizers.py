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


def _scalefactor(counts):
    """
    Calculate DESeq scaling factor per sample.
    """
    # masked array to discard inf, -inf, and nan
    ma = np.ma.masked_invalid(counts)
    return np.exp(np.ma.median(ma))


def norm_deseq(exp_obj):
    """
    Normalize by DESeq scaling factor, which is computed as the median of
    the ratio, for each row (gene), of its read count over its geometric
    mean across all samples. Return new counts dataframe.

    Details:
    --------
    http://genomebiology.com/2010/11/10/R106

    Parameters:
    -----------
    - exp_obj: experiment object. Normalized by DESeq scaling factor.
    """
    df = exp_obj.counts_df.copy()
    # log of counts
    lg = df.apply(np.log)
    # per sample: exponential(median(log(counts) - geometric mean))
    sf = lg.sub(lg.mean(axis=1), axis=0).apply(_scalefactor, axis=0)
    # apply scaling
    df = df.div(sf, axis=1)
    return df
