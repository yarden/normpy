##
## Representation of an RNA-Seq sample and its associated counts files
##

class Sample:
    def __init__(self, label, counts_fname):
        self.label = label
        self.counts_fname = counts_fname


