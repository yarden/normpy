##
## Representation of an RNA-Seq sample and its associated counts files
##
from collections import defaultdict

class Sample:
    """
    An RNA-Seq sample.
    """
    def __init__(self, label, counts_fname,
                 lib_size=None):
        self.label = label
        self.counts_fname = counts_fname
        self.lib_size = lib_size

    def get_lib_size():
        if self.lib_size is not None:
            return self.lib_size
        else:
            pass
        

