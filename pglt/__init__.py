#! /bin/usr/env python
# D.J. Bennett
# 26/05/2014

# POPULATE NAMESPACE
# http://stackoverflow.com/questions/4519127/setuptools-package-data-folder-location
import os
_ROOT = os.path.abspath(os.path.dirname(__file__))
_PARS = os.path.join(_ROOT, 'parameters.csv')
_GPARS = os.path.join(_ROOT, 'gene_parameters.csv')
import tools
import stages
# add stages -- a dictionary of stage functions -- to Stager
tools.system_tools.Stager.STAGES = stages.STAGES
# add default pars roots to Runner
tools.system_tools.Runner._pars = _PARS
tools.system_tools.Runner._gpars = _GPARS
# version (https://www.python.org/dev/peps/pep-0440/)
__version__ = open(os.path.join(os.path.dirname(__file__),
                                'VERSION.txt')).read().strip()
# docstring
__doc__ = open(os.path.join(os.path.dirname(__file__),
                            'DESCRIPTION.txt')).read()
del os
