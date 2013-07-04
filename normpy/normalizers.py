##
## Normalizers
##
## Methods:
##  UQ (Upper Quartile)
##  TMM (Trimmed Mean of M values)
##  TC (Total counts)
##  LOWESS (on Total counts)
##  [optional] RPKM (RPKM)
import numpy as np
import scipy
import pandas

import normpy
import normpy.utils as utils

try:
    # Import lowess function from statsmodels -- alternative to R
    from statsmodels.nonparametric.smoothers_lowess import lowess
except:
    raise Exception, "statsmodels not available?"

try:
    import rpy2
    from rpy2.robjects import r
    import rpy2.robjects as robj
    import rpy2.robjects.numpy2ri
    from rpy2.robjects.packages import importr
    rpy2.robjects.numpy2ri.activate()
    from rpy2.robjects.lib import grid
    from rpy2.robjects import r, Formula
    py2ri_orig = rpy2.robjects.conversion.py2ri
except:
    raise Exception, "Cannot import rpy2."



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


def run_lowess(X, Y,
               frac=0.75,
               missing="none"):
    """
    Y ~ X lowess.

    Parameters:
    -----------

    X: X values
    Y: Y values
    frac: fraction of data used to estimate each y-value.
    missing: how to handle missing values (by default "drop" them).
    """
    X[utils.where_null(X)] = np.nan
    Y[utils.where_null(Y)] = np.nan
    # Lowess takes Y values first
    fitted_Y = lowess(Y, X,
                      return_sorted=False,
                      frac=frac,
                      missing=missing)
    return fitted_Y
        

def R_norm_ma_lowess(exp_obj, pairs,
                     frac=0.75,
                     value_method="tc"):
    """
    Use R lowess.
    """
    # Total counts normalized dataframe of counts
    if value_method == "tc":
        df = norm_tc(exp_obj)
        df = utils.drop_null_rows(df, columns=utils.flatten(pairs))
    else:
        raise Exception, "Unknown value method %s" %(value_method)
    num_pairs = len(pairs)
    print "LOWESS MA normalizing %d pairs" %(num_pairs)
    norm_df = df.copy()
    norm_df = utils.drop_null_rows(norm_df, columns=utils.flatten(pairs))
    for pair in pairs:
        x, y = df[pair[0]], df[pair[1]]
        M = np.log2(x) - np.log2(y)
        # A = average intensity 1/2(XY)
        A = 0.5 * (np.log2(x) + np.log2(y))
        # Fit M ~ A
        corrected_M, corrected_factor = R_run_loess(A, M)
        corrected_x = 2**((2*A + corrected_M)/2.)
        corrected_y = 2**((2*A - corrected_M)/2.)
        del norm_df[pair[0]]
        del norm_df[pair[1]]
        norm_df[pair[0]] = corrected_x
        norm_df[pair[1]] = corrected_y
    return norm_df, df
    

def R_run_ma_loess(x, y):
    """
    Run MA-based loess normalization on X and Y. Computes

      M = log(X/Y)
      A = 0.5 * log(X*Y)

    Fits loess regression M ~ A and corrects X and Y accordingly.

    Assumes input X and Y values are non-logged.
    """
    M = np.log2(x) - np.log2(y)
    # A = average intensity 1/2(XY)
    A = 0.5 * (np.log2(x) + np.log2(y))
    # Fit M ~ A
    corrected_m, correction_factor = R_run_loess(A, M)
    corrected_x = 2**((2*A + corrected_m)/2.)
    corrected_y = 2**((2*A - corrected_m)/2.)
    return corrected_x, corrected_y


def R_run_loess(x, y, span=0.75):
    """
    Predict y as function of x. Takes two numpy vectors.
    """
    # Ensure that Inf/-Inf values are substituted
    x[utils.where_null(x)] = robj.NA_Real
    y[utils.where_null(x)] = robj.NA_Real
    data = robj.DataFrame({"x": x, "y": y})
    loess_fit = r.loess("y ~ x", data=data, span=span,
                        family="symmetric")
    correction_factor = np.array(list(r.predict(loess_fit, x)))
    corrected_y = \
        np.array(list(y)) - correction_factor
    return corrected_y, correction_factor

    
def norm_ma_lowess(exp_obj, pairs,
                   frac=0.75,
                   missing="none",
                   value_method="tc"):
    """
    Run MA-based loess normalization on total
    count normalized expression values. Computes

      M = log(X/Y)
      A = 0.5 * log(X*Y)

    where X, Y are total count normalized expression values
    for sample X and sample Y, respectively.

    Fits loess regression M ~ A and corrects X and Y accordingly.

    Assumes input X and Y values are non-logged.

    Parameters:
    -----------

    exp_obj: Experiment object
    pairs: Pairs to consider when doing the normalization
    frac: frac parameter to lowess() 
    missing: missing parameter to lowess()
    value_method: Method to use to get values for each sample.
    If 'tc', then normalize counts by total counts and then use
    that for lowess.

    Returns a normalized dataframe of values followed by the
    dataframe of values used in the normalization (the ones obtained
    by 'value_method', since lowess is not done on the counts.)
    """
    # Total counts normalized dataframe of counts
    if value_method == "tc":
        df = norm_tc(exp_obj)
        df = utils.drop_null_rows(df, columns=utils.flatten(pairs))
    else:
        raise Exception, "Unknown value method %s" %(value_method)
    num_pairs = len(pairs)
    print "LOWESS MA normalizing %d pairs" %(num_pairs)
    norm_df = df.copy()
    norm_df = utils.drop_null_rows(norm_df, columns=utils.flatten(pairs))
    for pair in pairs:
        x, y = df[pair[0]], df[pair[1]]
        M = np.log2(x) - np.log2(y)
        # A = average intensity 1/2(XY)
        A = 0.5 * (np.log2(x) + np.log2(y))
        # Fit M ~ A
        corrected_M = run_lowess(A, M,
                                 frac=frac,
                                 missing=missing)
        corrected_x = 2**((2*A + corrected_M)/2.)
        corrected_y = 2**((2*A - corrected_M)/2.)
        del norm_df[pair[0]]
        del norm_df[pair[1]]
        norm_df[pair[0]] = corrected_x
        norm_df[pair[1]] = corrected_y
    return norm_df, df


def norm_rpkm():
    """
    RPKM [optional].
    """
    pass


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
