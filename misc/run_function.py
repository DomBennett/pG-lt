#! /bin/usr/env python
# D.J. Bennett
# 20/11/2014
# The run function:
#  Takes a setup folder containing a list of names and sequence data, generates
#  a distribution of trees.
#  Requires installation of pglt package.

# PACKAGES
import os
from pglt.tools.system_tools import Stager


# RUN
def run(folder):
    # get working directory
    wd = os.path.join(os.getcwd(), folder)
    # run stages 3 (alignment) and 4 (phylogeny)
    Stager.run_all(wd=wd, stages=['3', '4'])
