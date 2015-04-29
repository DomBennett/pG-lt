#! /bin/usr/env python
# D.J. Bennett
# 06/11/2014
"""
pglt special tools
"""

# PACKAGES
import re
import os
import time
import sys
import subprocess
import shutil
import pickle
from tabulate import tabulate


# CLASSES
class Reseter(object):
    '''Reseter class : change settings in pG-lt folders'''
    options_msg = '''
    Options:
    1 - Reset to previous stage for all in `folders`
    2 - Change general parameters for all in `folders`
    3 - Change gene parameters for all in `folders`
    4 - Specify `folders` and return to these options

    By default, `folders` is a list of all pG-lt generated folders in current
    directory. To specify the folders for which you want to change settings,
    select option 4.

    Current directory: {0}

    To exit, press ctrl+c at any time.
    '''.format(os.getcwd())
    parameter_keys = ['nseqs', 'naligns', 'nphylos', 'thoroughness',
                      'parentid', 'outgroupid', 'maxtrys', 'rttstat',
                      'threads']
    # numerical only at this stage
    gene_parameter_keys = ['minlen', 'maxlen', 'mingaps', 'minoverlap',
                           'maxtrys', 'minseedsize', 'maxseedsize',
                           'maxseedtrys']
    gene_keys = ['NADHs1', 'NADHs2', 'NADHs3', 'NADHs4', 'NADHs5', 'NADHs6',
                 'NADHs7', 'NADHs8', 'NADHs9', 'NADHs10', 'NADHs11', 'NADHs12',
                 'COI', 'CYTB', 'COII', '12S', '16S', '18S', '28S', 'rbcl',
                 'matk']

    def __init__(self):
        folders = os.listdir(os.getcwd())
        avoid_pattern = "^\.|^log\.txt$|resolved_names|README|tempfiles"
        folders = [e for e in folders if not re.search(avoid_pattern, e)]
        self.folders = folders
        self.backup_folders = folders

    def _optionPrint(self, options):
        for option in options:
            print('    {0}'.format(option))

    def _deleteStageFolder(self, study_folder, stage_folder):
        '''Delete stage folder'''
        stage_folder = os.path.join(study_folder, stage_folder)
        if os.path.isdir(stage_folder):
            shutil.rmtree(stage_folder)

    def _readPickledFile(self, folder, filename):
        '''Read in pickled tempfile'''
        filepath = os.path.join(folder, 'tempfiles', filename)
        if not os.path.isfile(filepath):
            return None
        with open(filepath, "rb") as file:
            pickled = pickle.load(file)
        return(pickled)

    def _writePickledFile(self, folder, filename, pickled):
        '''Write out pickled tempfile'''
        filepath = os.path.join(folder, 'tempfiles', filename)
        with open(filepath, "wb") as file:
            pickle.dump(pickled, file)

    def _resetstage(self):
        """Reset `folders` to previous stage"""
        print('-'*70)
        while True:
            print('Reset `folders` to previous stages.')
            stage = raw_input('Enter stage (1-3): ')
            if int(stage) <= 3 and int(stage) > 0:
                break
            else:
                print('Invalid stage!')
        counter = 0
        stagenames = {'1': '1_names', '2': '2_download', '3': '3_alignment',
                      '4': '4_phylogeny'}
        stages = range(int(stage), 5)
        for folder in self.folders:
            # if progress exists, read in progress.p, delete stage folders
            progress = self._readPickledFile(folder=folder,
                                             filename='progress.p')
            if progress:
                for s in stages:
                    if progress[str(s)] != 'not run':
                        progress[str(s)] = 'not run'
                        self._deleteStageFolder(folder, stagenames[str(s)])
                self._writePickledFile(folder=folder, filename='progress.p',
                                       pickled=progress)
                counter += 1
        # reset global progress
        progress = self._readPickledFile(folder='', filename='progress.p')
        for s in stages:
            progress[str(s)] = 'not run'
        self._writePickledFile(folder='', filename='progress.p',
                               pickled=progress)
        print('    [{0}] folders reset to stage [{1}]'.format(counter, stage))

    def _resetparameters(self):
        '''Change key's value in all paradicts in `folders`'''
        print('-'*70)
        while True:
            print('Change parameters for all `folders`. Available options:')
            self._optionPrint(self.parameter_keys)
            key = raw_input('Enter parameter name of setting to change: ')
            if key not in self.parameter_keys:
                print('Invalid parameter name!')
            else:
                break
        value = raw_input('Enter new parameter setting: ')
        counter = 0
        if key == 'threads':
            for folder in self.folders:
                threads = self._readPickledFile(folder=folder,
                                                filename='threads.p')
                if threads:
                    threads = value
                    self._writePickledFile(folder=folder, filename='threads.p',
                                           pickled=threads)
                    counter += 1
        else:
            for folder in self.folders:
                paradict = self._readPickledFile(folder=folder,
                                                 filename='paradict.p')
                if paradict:
                    paradict[key] = value
                    self._writePickledFile(folder=folder,
                                           filename='paradict.p',
                                           pickled=paradict)
                    counter += 1
        print('    [{0}] set to [{1}] for [{2}] folders'.
              format(key, value, counter))

    def _resetgeneparameters(self):
        '''Change key's value in a genedict in `folders`'''
        print('-'*70)
        while True:
            print('Enter name of gene to change. Available options:')
            self._optionPrint(self.gene_keys)
            gene = raw_input('Enter name: ')
            if gene in self.gene_keys:
                break
            else:
                print('Invalid gene name!')
        while True:
            print('Enter parameter name of setting to change. Available options:')
            self._optionPrint(self.gene_parameter_keys)
            key = raw_input('Enter name: ')
            if key in self.gene_parameter_keys:
                break
            else:
                print('Invalid parameter name!')
        value = raw_input('Enter new parameter setting: ')
        counter = 0
        for folder in self.folders:
            genedict = self._readPickledFile(folder=folder,
                                             filename='genedict.p')
            if genedict:
                genedict[gene][key] = value
                self._writePickledFile(folder=folder, filename='genedict.p',
                                       pickled=genedict)
                counter += 1
        print('    [{0}] set to [{1}] for [{2}] gene for [{3}] folders'.
              format(key, value, gene, counter))

    def _setfolders(self):
        '''Change `folders`'''
        print('-'*70)
        user_folders = raw_input('Enter folder names separated by spaces: ')
        user_folders = user_folders.split(' ')
        for folder in user_folders:
            if folder not in self.backup_folders:
                print('[{0}] not found -- did you spell it correctly?'.
                      format(folder))
        self.folders = user_folders

    def run(self):
        '''Run reset'''
        try:
            while True:
                print('-'*70)
                print('\n{0} RESET MODE {0}'.format(' '*29))
                print('-'*70)
                print(self.options_msg)
                option = str(raw_input('Enter option number (1-4): '))
                if '1' == option:
                    self._resetstage()
                if '2' == option:
                    self._resetparameters()
                if '3' == option:
                    self._resetgeneparameters()
                if '4' == option:
                    self._setfolders()
                # TODO: reset the ones that failed
        except KeyboardInterrupt:
            sys.exit('\nExiting reset mode ....')


# FUNCTIONS
def timeit(func, **kwargs):
    """Time a function (platform independent)"""
    # http://stackoverflow.com/questions/10884751/python-timeit-and-program-output
    if sys.platform == "win32":
        # On Windows, the best timer is time.clock()
        default_timer = time.clock
    else:
        # On most other platforms the best timer is time.time()
        default_timer = time.time
    t0 = default_timer()
    output = func(**kwargs)
    t1 = default_timer()
    return output, t1-t0


def getThreads(wd=None):
    """Find number of cores on machine (platform independent). If wd, \
search wd for pickled threads.p."""
    if wd:
        filepath = os.path.join(wd, 'threads.p')
        if os.path.isfile(filepath):
            with open(filepath, "rb") as file:
                nthreads = pickle.load(file)
        else:
            nthreads = 1
        return nthreads
    if sys.platform == 'darwin':
        # find threads for a mac
        cmd = "system_profiler SPHardwareDataType | grep Cores | awk \
'{print $5}'"
    elif sys.platform == 'linux2':
        # find threads for a unix
        cmd = "cat /proc/cpuinfo | grep processor | tail -1 | awk '{print $3}'"
    else:
        # TODO: work out n cores for windows
        # TODO: add a threads parameter to parameters.csv
        return 1
    s = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                         shell=True)
    nthreads = s.stdout.readline()
    try:
        nthreads = int(nthreads)
    except:
        nthreads = 1
    return nthreads


def clean():
    """Remove all pglt files and folders"""
    # remove log.txt in parent folder
    try:
        os.remove('log.txt')
    except OSError:
        pass

    # remove tempfiles
    try:
        shutil.rmtree('tempfiles')
    except OSError:
        pass

    # remove resolved_names
    try:
        shutil.rmtree('resolved_names')
    except OSError:
        pass

    # go through all subdirs and remove files in list
    files_to_remove = ['info.txt', 'log.txt']
    folders = os.listdir(os.getcwd())
    while folders:
        temp_folder = folders.pop()
        temp_files_to_remove = files_to_remove[:]
        while temp_files_to_remove:
            temp_file_to_remove = temp_files_to_remove.pop()
            try:
                os.remove(os.path.join(temp_folder, temp_file_to_remove))
            except OSError:
                pass

    # go through all subdirs and remove folders in list
    folders_to_remove = ['1_names', '2_download', '3_alignment', '4_phylogeny',
                         'tempfiles']
    folders = os.listdir(os.getcwd())
    while folders:
        temp_folder = folders.pop()
        temp_folders_to_remove = folders_to_remove[:]
        while temp_folders_to_remove:
            temp_folder_to_remove = temp_folders_to_remove.pop()
            try:
                shutil.rmtree(os.path.join(temp_folder, temp_folder_to_remove))
            except OSError:
                pass


def check(stage, stagecounter, directory):
        '''Check stage'''
        progress_path = os.path.join(directory, 'tempfiles', 'progress.p')
        if not os.path.isfile(progress_path):
            stagecounter[stage]['not run'] += 1
        else:
            with open(progress_path, "rb") as file:
                progress = pickle.load(file)
            stagecounter[stage][progress[stage]] += 1
        return stagecounter


def stats():
    """Print success status of folders"""
    # counter dicts
    statecounter = {'not run': 0, 'failed': 0, 'success': 0}
    stagecounter = {'1': statecounter.copy(), '2': statecounter.copy(),
                    '3': statecounter.copy(), '4': statecounter.copy()}
    consensus_counter = 0
    # get all folders
    folders = os.listdir(os.getcwd())
    avoid_pattern = "^\.|^log\.txt$|resolved_names|README|tempfiles"
    folders = [e for e in folders if not re.search(avoid_pattern, e)]
    stages = ['1', '2', '3', '4']
    for folder in folders:
        for stage in stages:
            stagecounter = check(stage=stage, stagecounter=stagecounter,
                                 directory=os.path.join(os.getcwd(), folder))
        if os.path.isfile(os.path.join(os.getcwd(), folder, '4_phylogeny',
                                       'consensus.tre')):
            consensus_counter += 1
    print ('\n{0}/{1} folders with phylogenetic trees....\n'.
           format(consensus_counter, len(folders)))
    headers = ['Stage', 'Not run', 'Failed', 'Success']
    rows = []
    for stage in stages:
        sc = stagecounter[stage]
        rows.append([stage, sc['not run'], sc['failed'], sc['success']])
    print(tabulate(rows, headers, tablefmt="simple"))
