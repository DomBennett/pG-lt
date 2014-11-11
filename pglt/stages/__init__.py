#! /usr/bin/env python
# D.J. Bennett
# 26/04/2014
"""
pglt stages
"""
# populate STAGES dictionary with run stage functions
#  and output dirs
import names_stage
import download_stage
import alignment_stage
import phylogeny_stage
STAGES = {'1': (names_stage.run, '1_names'),
          '2': (download_stage.run, '2_download'),
          '3': (alignment_stage.run, '3_alignment'),
          '4': (phylogeny_stage.run, '4_phylogeny')}
