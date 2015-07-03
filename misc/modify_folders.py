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
    gdpath = os.path.join(folder, 'tempfiles', 'genedict.p')
    if os.path.isfile(gdpath):
        counter += 1
        with open(gdpath, "rb") as file:
            genedict = pickle.load(file)
        rna_genes = ['12S', '16S', '18S', '28S']
        genes = [e for e in genedict.keys() if e not in rna_genes]
        for g in genes:
            genedict[g]['maxgaps'] = 10
        for rg in rna_genes:
            genedict[rg]['maxgaps'] = 20
        with open(gdpath, "wb") as file:
            pickle.dump(genedict, file)
print '[{0}] folders'.format(counter)
