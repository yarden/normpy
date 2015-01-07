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
    def __init__(self, counts, samples,
                 gene_id_col="gene_id",
                 lib_sizes=None,
                 reads_type="single",
                 exp_type="vanilla",
                 extra_cols=[]):
        """
        Setup the experiments.

        Parameters:
        ----------
        - counts: if string, then treated as a filename with counts for each sample.
        if of type pandas.DataFrame, then used as dataframe of counts.
        - samples: samples in the experiment. Mapping from sample labels
        to columns in a counts filename.
        - gene_id_col: column that identifies the gene id. Must be unique
        in the counts dataframe.
        - lib_sizes: library sizes. Mapping from sample column to size of library.
        If not given, this is computed as the sum of counts in each sample. 
        - reads_type: single-end or paired-end. 
        - exp_type: The experimental design type. (Not used
        for now.)
        - extra_cols: optional list of extra columns to import from
        data frame into normalized column
        """
        self.counts = counts
        self.samples = samples
        self.num_samples = len(self.samples)
        self.gene_id_col = gene_id_col
        self.lib_sizes = lib_sizes
        self.reads_type = reads_type
        self.exp_type = exp_type
        self.extra_cols = extra_cols
        # dataframe with counts for all samples
        self.counts_df = None
        self.load_counts()
        

    def load_counts(self, sep="\t"):
        """
        Load counts for all samples into pandas dataframe.
        """
        if type(self.counts) == str:
            if not os.path.isfile(self.counts):
                raise Exception, "Cannot find counts file %s" %(self.counts)
            counts_fname = self.counts
            self.counts_df = pandas.read_table(counts_fname, sep=sep)
        else:
            # Assume we were given a dataframe
            self.counts_df = self.counts
        # Set gene id column to be the index unless it is
        # not found in the table, in which case assume
        # it is already indexed by that column
        if self.gene_id_col in self.counts_df:
            self.counts_df = self.counts_df.set_index(self.gene_id_col)
        # Select only the columns relevant to the samples, plus
        # the gene id column
        relevant_cols = self.extra_cols + self.samples.values()
        self.counts_df = self.counts_df[relevant_cols]


    def __repr__(self):
        return self.__str__()


    def __str__(self):
        str_obj = \
            "Experiment(type=%s, num_samples=%d, reads_type=%s)" \
            %(self.exp_type,
              self.num_samples,
              self.reads_type)
        return str_obj

