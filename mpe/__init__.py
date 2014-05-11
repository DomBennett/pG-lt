# http://stackoverflow.com/questions/4519127/setuptools-package-data-folder-location
import os

_ROOT = os.path.abspath(os.path.dirname(__file__))
_PARS = os.path.join(_ROOT,'parameters.csv')
_GPARS = os.path.join(_ROOT,'gene_parameters.csv')
from mpe import tools
from mpe import stages

del os