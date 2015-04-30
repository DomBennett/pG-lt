#! /bin/usr/env python
# D.J. Bennett
# 26/05/2014

# POPULATE NAMESPACE
# http://stackoverflow.com/questions/4519127/setuptools-package-data-folder-location
import os
import pickle
_ROOT = os.path.abspath(os.path.dirname(__file__))
_PARS = os.path.join(_ROOT, 'parameters.csv')
_GPARS = os.path.join(_ROOT, 'gene_parameters.csv')
with open(os.path.join(_ROOT, 'dependencies.p'), "rb") as file:
    depsdict = pickle.load(file)
_RAXML = depsdict['raxml']
_MAFFT = depsdict['mafft']
_MAFFTQ = depsdict['mafftq']
_MAFFTX = depsdict['mafftx']
_BLASTN = depsdict['blastn']
del depsdict
import tools
import stages
# add stages -- a dictionary of stage functions -- to Stager
tools.system_tools.Stager.STAGES = stages.STAGES
# add default pars roots to Runner
tools.system_tools.Runner._pars = _PARS
tools.system_tools.Runner._gpars = _GPARS
# version (https://www.python.org/dev/peps/pep-0440/)
__version__ = 'v0.1-beta'
__year__ = '2015'
# docstring
__doc__ = '''
pG-lt (or phyloGenerator-light) is a pipeline for the automated generation of
phylogenies through `Mass Phylogeny Estimation`. It is adapted from
phyloGenerator (C) 2013 and was written by D.J. Bennett with additional help
and suggestions from W.D. Pearse, L. Hudson and A. Purvis.

It makes use of external programs for phylogeny generation and bioinformatics,
these are: RAxML (Copyright (C) Stamatakis 2013), MAFFT (Copyright (C) 2013
Kazutaka Katoh) NCBI's standalone BLAST suite 2.2.29+ and online API services
(Copyright NCBI (C) 2009). It also uses the following python packages:
Biopython (Copyright Cook (C) 2009), Dendropy (Copyright Sukumaran and Holder
(C) 2010) and Taxon Names Resovler (Copyright (C) Bennett 2014).

For details on how to use pG-lt, please refer to its wiki:
`https://github.com/DomBennett/pG-lt/wiki`
For any questions or comments, feel free to email:
`dominic.john.bennett@gmail.com`.

Copyright (C) {0} Dominic John Bennett

This program is free software; you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation; either version 2 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
Street, Fifth Floor, Boston, MA 02110-1301 USA.
'''.format(__year__)
del os
del pickle
