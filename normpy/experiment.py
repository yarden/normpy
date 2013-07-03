##
## Representation of an Experiment and a set of
## samples.  Keep it simple.
##
from collections import OrderedDict

class Experiment:
    """
    Representation of an Experiment, which involves a set
    of Samples and technical parameters of the experiment
    (experimental design, the type of sequencing reads, etc.)
    """
    def __init__(self, samples,
                 reads_type="single",
                 exp_type=None):
        """
        Initialize experiments.

        - samples: samples in the experiment
        - reads_type: single-end or paired-end
        - exp_type: The experimental design type. Not used
        for now.
        """
        self.samples = samples
        self.reads_type = reads_type
        self.exp_type = exp_type


