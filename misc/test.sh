#!/bin/bash
# Install, setup, test, send results to coveralls.io

# clone
git clone https://github.com/DomBennett/pG-lt.git
# install
cd pglt
python setup.py install
# set deps in local copy
pglt_set_dependencies.py -mafft $(which mafft) -raxml $(which raxml) -blast $(which blastn) -mafftq $(which mafft-qinsi) -mafftx $(which mafft-xinsi) --local --overwrite
# run tests and coverage
coverage run setup.py test
# report
coverage report
coveralls
