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


def _sf_deseq(counts):
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
    sf = lg.sub(lg.mean(axis=1), axis=0).apply(_sf_deseq, axis=0)
    # apply scaling
    df = df.div(sf, axis=1)
    return df

def _sf_tmm(obs, ref, log_ratio_trim=0.3, sum_trim=0.05, weighting=True, a_cutoff=-1e10):
    """
    Called by `norm_tmm`.
    """
    if all(obs == ref):
        return 1
    
    obs_sum = obs.sum()
    ref_sum = ref.sum()

    # log ration of expression accounting for library size
    lr = np.log2((obs / obs_sum) / (ref / ref_sum))
    # absolute expression
    ae = (np.log2(obs / obs_sum) + np.log2(ref / ref_sum)) / 2
    # estimated asymptotic variance
    v = (obs_sum - obs) / obs_sum / obs + (ref_sum - ref) / ref_sum / ref
    # create mask
    m = np.isfinite(lr) & np.isfinite(ae) & (ae > -1e10)
    # drop the masked values
    lr = lr[m]
    ae = ae[m]
    v = v[m]
    assert len(lr) == len(ae) == len(v)
    
    n = len(lr)
    lo_l = np.floor(n * log_ratio_trim) + 1
    hi_l = n + 1 - lo_l
    lo_s = np.floor(n * sum_trim) + 1
    hi_s = n + 1 - lo_s
    k = ((lr.rank(method="first") >= lo_l) & (lr.rank(method="first") <= hi_l)) \
            & ((ae.rank(method="first") >= lo_s) & (ae.rank(method="first") <= hi_s))

    if weighting:
        return 2**(sum(lr[k] / v[k]) / sum(1 / v[k]))
    else:
        return 2**(lr[k].mean())

def _sf_q(df, q=0.75):
    """
    Parameters:
    -----------
    - df: zeroed rows removed
    - q: quartile
    """
    lib_size = df.sum()
    y = df.T.div(lib_size, axis=0).T
    # fill nans with 0
    y = y.dropna(how="all").fillna(0)
    y = y.quantile(q)
    # factors multiple to one
    sf = y.div(np.exp(np.mean(np.log(y))))
    return sf

def norm_tmm(exp_obj, ref_col=None, log_ratio_trim=0.3, sum_trim=0.05, weighting=True, a_cutoff=-1e10):
    """
    Trimmed Mean of M-values (TMM) is the weighted mean of log ratios between
    this test and the reference, after exclusion of the most expressed genes
    and the genes with the largest log ratios.
    
    Parameters:
    -----------
    - exp_obj: experiment object.
    - ref_col: reference column from which to scale others.
    - log_ratio_trim: amount of trim to use on log-ratios.
    - sum_trim: amount of trim to use on combined absolute values.
    - weighting: whether to compute weights.
    - a_cutoff: cutoff of absolute expression values.
    """
    df = exp_obj.counts_df.copy()
    # remove zeros
    nz = df.where(df > 0)
    nz = nz.dropna(how="all").fillna(0)
    # reference column
    if ref_col is None:
        # quantile factors
        sf_q = _sf_q(nz)
        ref_col = (abs(sf_q - np.mean(sf_q))).idxmin()
    # try:
    kwargs = {"ref":nz[ref_col],
                "log_ratio_trim":log_ratio_trim,
                "sum_trim":sum_trim,
                "weighting":weighting,
                "a_cutoff":a_cutoff}
    # except KeyError:
        # revert back to auto?
    sf_tmm = nz.apply(_sf_tmm, **kwargs)
    # apply scaling
    df = df.div(sf_tmm, axis=1)
    return df

def norm_q(exp_obj, q=0.75):
    """
    Ported from edgeR and still needs to be validated. Also, maybe compare
    edgeR method to limma implementation.
    
    Quantile normalization.
    
    Parameters:
    -----------
    - exp_obj: experiment object.
    - q: quantile.
    """
    df = exp_obj.counts_df.copy()
    # remove zeros
    nz = df.where(df > 0)
    nz = nz.dropna(how="all").fillna(0)
    sf_q = _sf_q(nz, q)
    # apply scaling
    df = df.div(sf_q, axis=1)
    return df
