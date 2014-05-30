#! /bin/usr/env python
## D.J. Bennett
## 26/05/2014
"""
MPE is a pipeline for the automated generation of phylogenies through 'Mass
Phylogeny Estimation'. This program is built on top of phyloGenerator (C) 2013
and was written by D.J. Bennett with additional help from W.D. Pearse and L. Hudson.
This program makes use of external programs for phylogeny generation and bioinformatics
these are: RAxML (Copyright (C) Stamatakis 2013) , MAFFT (Copyright (C) 2013 Kazutaka
Katoh) the NCBI's standalone BLAST suite 2.2.29+ and online API services
 (Copyright NCBI (C) 2009). It also uses a variety of python packages including:
 Biopython (Copyright Cook (C) 2009) and Dendropy (Copyright Sukumaran and Holder (C)
 2010).

Copyright (C) 2014  Dominic John Bennett

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
"""

## Populate namespace
# http://stackoverflow.com/questions/4519127/setuptools-package-data-folder-location
import os
_ROOT = os.path.abspath(os.path.dirname(__file__))
_PARS = os.path.join(_ROOT,'parameters.csv')
_GPARS = os.path.join(_ROOT,'gene_parameters.csv')
import tools
import stages
# add stages -- a dictionary of stage functions -- to Stager
tools.system.Stager.STAGES = stages.STAGES
del os