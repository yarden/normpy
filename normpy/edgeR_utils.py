##
## Utilities for calling edgeR
##
import os
import sys
import time

import pandas

from collections import OrderedDict

import normpy
import normpy.utils as utils

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


class EdgeR:
    def __init__(self, exp_obj):
        """
        Parameters:
        -----------
        - exp_obj: an Experiments object to work with
        """
        load_edgeR()
        self.exp_obj = exp_obj
        

    def norm_expr_vals(self, ref_col, method="TMM"):
        """
        Normalize expression values relative to a
        reference sample (using TMM normalization).

        Parameters:
        -----------
        - ref_col: Reference column
        - method: Method to use for normalization, e.g. 'TMM'
        """
        # Take only numeric values from dataframe
        r_counts = conversion_pydataframe(self.exp_obj.counts_df)
        r_dge = r.DGEList(r_counts)
        # Calculate normalization factors
        r_dge = r.calcNormFactors(r_dge,
                                  refColumn=ref_col,
                                  method=method)
        print r_dge[0:10]
        # Get counts per million
        r_cpm_result = r.cpm(r_dge)
        print "Counts per million: "
        print r_cpm_result[0:10]
        return r_dge, r_cpm_result
        

def load_edgeR():
    """
    Load edgeR library into R.
    """
    robj.r("library(edgeR)")


def conversion_pydataframe(obj):
    """
    Convert pandas DataFrame or python object to an R dataframe/object.
    """
    if isinstance(obj, pandas.core.frame.DataFrame):
        od = OrderedDict()
        for name, values in obj.iteritems():
            if values.dtype.kind == 'O':
                od[name] = rpy2.robjects.vectors.StrVector(values)
            else:
                od[name] = rpy2.robjects.conversion.py2ri(values)
        return rpy2.robjects.vectors.DataFrame(od)
    else:
        return py2ri_orig(obj)


