#! /usr/bin/env python
"""
Run me in a pG-lt setup folder to run individual stages in a python session
"""

# LIBS
import logging
import pglt

# SETUP LOGGER
logger = logging.getLogger('')  # this creates a root logger
logger.setLevel(logging.DEBUG)  # print debug messages and above
console = logging.StreamHandler()  # create console stream
console.setFormatter(logging.Formatter('%(message)s'))  # simple message
logger.addHandler(console)  # add to logger

# RUN (names, download, alignment of phylogeny)
pglt.stages.alignment_stage.run()
