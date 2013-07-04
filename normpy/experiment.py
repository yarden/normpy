##
## Representation of an Experiment and a set of
## samples.  Keep it simple.
##
import os
import sys
import time

import numpy as np
import pandas

from collections import OrderedDict


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


class Experiment:
    """
    Representation of an Experiment, which involves a set
    of Samples and technical parameters of the experiment
    (experimental design, the type of sequencing reads, etc.)
    """
    def __init__(self, counts_fname, samples,
                 gene_id_col="gene_id",
                 lib_sizes=None,
                 reads_type="single",
                 exp_type=None):
        """
        Setup the experiments.

        Parameters:
        ----------
        - counts_fname: filename with counts for each sample in the experiment
        - samples: samples in the experiment. Mapping from sample labels
        to columns in a counts filename.
        - gene_id_col: column that identifies the gene id. Must be unique
        in the counts dataframe.
        - lib_sizes: library sizes. Mapping from sample column to size of library.
        If not given, this is computed as the sum of counts in each sample. 
        - reads_type: single-end or paired-end. 
        - exp_type: The experimental design type. (Not used
        for now.)
        """
        self.counts_fname = counts_fname
        self.samples = samples
        self.gene_id_col = gene_id_col
        self.lib_sizes = lib_sizes
        self.reads_type = reads_type
        self.exp_type = exp_type
        # dataframe with counts for all samples
        self.counts_df = None
        self.load_counts()
        

    def load_counts(self, sep="\t"):
        """
        Load counts for all samples into pandas dataframe.
        """
        if not os.path.isfile(self.counts_fname):
            raise Exception, "Cannot find counts file %s" %(self.counts_fname)
        self.counts_df = pandas.read_table(self.counts_fname, sep=sep)
        # Select only the columns relevant to the samples, plus
        # the gene id column
        relevant_cols = [self.gene_id_col] + self.samples.values()
        self.counts_df = self.counts_df[relevant_cols]
        # Set gene id column to be the index
        self.counts_df = self.counts_df.set_index(self.gene_id_col)
