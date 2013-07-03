##
## Representation of an RNA-Seq sample and its associated counts files
##
from collections import defaultdict

class Sample:
    """
    An RNA-Seq sample.
    """
    def __init__(self, label, counts_fname):
        self.label = label
        self.counts_fname = counts_fname

