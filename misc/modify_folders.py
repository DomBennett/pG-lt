#! /usr/bin/env python
# D.J. Bennett
# 04/07/2015
"""
Loop through folders making changes
"""

# PACKAGES
import os
import pickle

# Here looping throgh, making changes to genedict
folders = os.listdir('.')
counter = 0
for folder in folders:
    parpath = os.path.join(folder, 'tempfiles', 'paradict.p')
    if os.path.isfile(parpath):
        counter += 1
        with open(parpath, "rb") as file:
            pardict = pickle.load(file)
        pardict['minspecies'] = '5'
        with open(parpath, "wb") as file:
            pickle.dump(pardict, file)
print '[{0}] folders'.format(counter)

# Here converting all progress dicts
folders = os.listdir('.')
counter = 0
for folder in folders:
    prgpath = os.path.join(folder, 'tempfiles', 'progress.p')
    if os.path.isfile(prgpath):
        counter += 1
        with open(prgpath, "rb") as file:
            progress = pickle.load(file)
        progress['1'] = 'not run'
        with open(prgpath, "wb") as file:
            pickle.dump(progress, file)
print '[{0}] folders'.format(counter)
