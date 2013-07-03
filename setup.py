from distutils.core import setup, Extension
import glob
import os
import sys

# This forces distutils to place the data files
# in the directory where the Py packages are installed
# (usually 'site-packages'). This is unfortunately
# required since there's no good way to retrieve
# data_files= from setup() in a platform independent
# way.
from distutils.command.install import INSTALL_SCHEMES
for scheme in INSTALL_SCHEMES.values():
        scheme['data'] = scheme['purelib']
        
setup(name = 'normpy',
      version = '0.1',
      description = 'Normalization methods for RNA-Seq data.',
      author = 'Yarden Katz,Rory Kirchner,Ryan Dale,Brad Chapman',
      author_email = 'yarden@mit.edu',
      packages = ['normpy'],
      package_dir={'normpy': 'normpy'},
      # Required modules
      install_requires = [
#          "matplotlib >= 1.1.0",
#          "matplotlib",
          "numpy >= 1.5.0",
          "scipy",
          "statsmodels",
          "pandas >= 0.8.1",
          ],
      platforms = 'ALL',
      keywords = ['bioinformatics', 'sequence analysis',
                  'gene expression', 'RNA-Seq'],
      classifiers = [
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
        ]
      )

