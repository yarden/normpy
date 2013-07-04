##
## Misc utilities for normpy
##
import os
import numpy as np
import pandas


def load_testdata(test_name):
    testdir = os.path.join("test", "data")
    if test_name == "pasilla":
        # Return pasilla counts filename
        pasilla_fname = \
            os.path.join(os.path.abspath(os.path.join(__file__, "..")),
                         testdir, "pasilla_gene_counts.tsv")
        if not os.path.isfile(pasilla_fname):
            raise Exception, "Cannot find pasilla file %s" %(pasilla_fname)
        return pasilla_fname
    else:
        raise Exception, "Unknown test name %s" %(test_name)


def get_plots_dir():
    plots_dir = \
        os.path.join(os.path.abspath(os.path.join(__file__, "..")),
                     "plots")
    return plots_dir


def where_na_like(l):
    """
    Return indices where array is NA-like
    """
    bool_index = np.array(map(lambda x: np.isinf(x) or \
                              pandas.isnull(x), l))
    return np.where(bool_index)[0]


def flatten(l):
    """
    Flatten a list.
    """
    return [item for sublist in l for item in sublist]


def where_null(l):
    """
    Remove null values including inf.
    """
    return l.apply(lambda x: x in [np.inf, -np.inf, np.nan])


def drop_null_rows(df, columns=None):
    new_df = select_df_rows(df,
                            lambda x: \
                            not x in [np.inf, -np.inf, np.nan],
                            columns=columns)[0]
    return new_df


def select_df_rows(df, cond,
                   columns=None,
                   how="any"):
    """
    Select rows from DataFrame where a condition
    holds.  cond is a lambda.
    """
    if columns is None:
        columns = df.columns
    # Apply condition result and get array as result
    cond_result = df[columns].applymap(cond).values
    result = None
    if how == 'any':
        result_ind = cond_result.any(1)
    elif how == 'all':
        result_ind = cond_result.all(1)
    elif type(how) == int:
        # apply to N many columns
        result_ind = cond_result.sum(axis=1) >= how
    else:
        raise Exception, "Do not know %s" %(how)
    result = df[result_ind]
    return result, result_ind


def make_dir(dirpath):
    if os.path.isfile(dirpath):
        print "Error: %s is a file!" %(dirpath)
        sys.exit(1)
    # Try to make the directory
    try:
        os.makedirs(dirpath)
    except OSError:
        pass

