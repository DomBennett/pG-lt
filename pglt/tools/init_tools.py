#! /bin/usr/env python
# D.J. Bennett
# 06/11/2014
"""
pglt init tools
"""

# PACKAGES
import os
import pickle
import csv
from system_tools import PrimingError
from special_tools import getThreads

# MESSAGES
ioerror_msg = "[{0}] file could not be opened in [{1}]. Check that \
it is not opened by another program"


# FUNCTIONS
def recordPars(paradict):
    """Return mpe parameters string"""
    record = '\nUsing the following parameters:\n'
    for key in paradict.keys():
        record += '    [{0}] = [{1}]\n'.format(key, paradict[key])
    return record


def recordGpars(genedict):
    """Return gene parameters string"""
    record = '\nUsing the following genes and gene parameters:\n'
    for gene in genedict.keys():
        record += '  Gene: [{0}]\n'.format(gene)
        for par in genedict[gene]:
            record += '    [{0}] = [{1}]\n'.format(par, genedict[gene][par])
    return record


def prime(directory, arguments):
    """Write pickle files, print arguments"""
    # Write pickle files
    with open(os.path.join(directory, ".genedict.p"), "wb") as file:
        pickle.dump(arguments['genedict'], file)
    with open(os.path.join(directory, ".paradict.p"), "wb") as file:
        pickle.dump(arguments['paradict'], file)
    with open(os.path.join(directory, ".terms.p"), "wb") as file:
        pickle.dump(arguments['terms'], file)
    # Print arguments and parameters to file
    record = 'Working with [{0}] names\n'.format(len(arguments['terms']))
    record += recordPars(arguments['paradict'])
    record += recordGpars(arguments['genedict'])
    with open(os.path.join(directory, 'info.txt'), 'wd') as file:
        file.write(record)
    # TODO return the starting stage
    # TODO create other functions for creating hidden files


def readInNames(directory):
    """Read names from text file in dir"""
    terms = []
    with open(os.path.join(directory, 'names.txt')) as names:
        for name in names:
            terms.append(name.strip())
    terms = [term for term in terms if not term == '']
    return terms


def readInGenePars(gpars_file, default_gpars_file):
    """Read gene_parameters.csv. Return list of dictionaries."""
    def _read(gpars_file, template, genes=None):
        # open csv file and replace parameters in template
        # if they are None. If genes specified, only read
        # rows for those genes.
        with open(gpars_file, 'rb') as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                if genes:
                    if not row['gene'] in genes:
                        continue
                temp = template.copy()
                for key in temp.keys():
                    if row[key]:
                        if temp[key] is None:
                            if key == 'names':
                                # for names, split into a list of syns
                                temp[key] = row[key].split(':')
                            else:
                                temp[key] = row[key]
                                genedict[row['gene']] = temp
        return genedict
    # check if file exists, else use default
    if not os.path.isfile(gpars_file):
        return readInGenePars(default_gpars_file, None)
    # genedicts
    genedict = {}
    # template of dict in genedict
    template = {'names': None, 'taxid': None, 'minlen': None, 'maxlen': None,
                'mingaps': None, 'minoverlap': None, 'minfails': None,
                'maxtrys': None, 'minseedsize': None, 'maxseedsize': None,
                'maxseedtrys': None, 'partition': None, 'type': None}
    # open file, read each row and fill in template
    genedict = _read(gpars_file, template)
    # if Nones, use defaults
    nones = False
    for gene in genedict.keys():
        for par in genedict[gene].keys():
            if genedict[gene][par] is None:
                nones = True
                break
    if nones:
        # run _read for defaults and limit to genes in genedict
        genedict = _read(default_gpars_file, template, genedict.keys())
    return genedict


def readInPars(pars_file, default_pars_file):
    """Read gene_parameters.csv. Return dictionary."""
    def _read(pars_file, paradict):
        # open csv, and replace all Nones
        with open(pars_file, 'rb') as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                if paradict[row["Parameter"]] is None:
                    paradict[row["Parameter"]] = row["Value"]
        return paradict
    # check if file exists, else use default
    if not os.path.isfile(pars_file):
        return readInPars(default_pars_file)
    # template
    paradict = {'nseqs': None, 'naligns': None, 'ntrees': None,
                'thoroughness': None, 'maxtrys': None, 'rttpvalue': None,
                'parentid': None, 'outgroupid': None, 'threads': None}
    # open file, read each row, extract value
    paradict = _read(pars_file, paradict)
    # if Nones remain, use default
    nones = False
    for key in paradict.keys():
        if paradict[key] is None:
            nones = True
            break
    if nones:
        paradict = _read(default_pars_file, paradict)
    # add threads if -1 is specified
    if int(paradict['threads']) == -1:
        paradict['threads'] = getThreads()
    return paradict


def sortArgs(directory, email, logger, default_pars_file, default_gpars_file):
    """Search for relevant files in dir, return list of arguments"""
    # find text file and read, raise error if fail
    try:
        terms = readInNames(directory)
    except IOError:
        logger.error(ioerror_msg.format('names.txt', directory))
        raise PrimingError()
    # find gene parameter file and read, raise error if fail
    try:
        genedict = readInGenePars(os.path.join(directory,
                                               'gene_parameters.csv'),
                                  default_gpars_file)
    except IOError:
        logger.error(ioerror_msg.format('gene_parameters.csv', directory))
        raise PrimingError()
    # find parameter file and read, raise error if fail
    try:
        paradict = readInPars(os.path.join(directory,
                                           'parameters.csv'),
                              default_pars_file)
    except IOError:
        logger.error(ioerror_msg.format('parameters.csv', directory))
        raise PrimingError()
    # add email to paradict
    paradict['email'] = email
    return {'terms': terms, 'genedict': genedict, 'paradict': paradict}
