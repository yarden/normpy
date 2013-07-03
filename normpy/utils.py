##
## Misc utilities for normpy
##
import os

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
