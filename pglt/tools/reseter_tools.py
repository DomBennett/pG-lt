#! /bin/usr/env python
# D.J. Bennett
# 21/05/2015
"""
pglt reseter tools
"""

# PACKAGES
import re
import os
import sys
import shutil
import pickle


# CLASSES
class Reseter(object):
    '''Reseter class : change settings in pG-lt folders'''
    verbose = True
    options_msg = '''
    Options:
    1 - Reset to previous stage for all in `folders`
    2 - Change general parameters for all in `folders`
    3 - Change gene parameters for all in `folders`
    4 - Specify `folders` and return to these options
    5 - Restore to defaults for all in `folders`

    By default, `folders` is a list of all pG-lt generated folders in current
    directory. To specify the folders for which you want to change settings,
    select option 4.

    Current directory: {0}

    To exit, press ctrl+c at any time.
    '''

    def __init__(self, paradict, genedict, wd=os.getcwd()):
        self.wd = wd
        folders = os.listdir(wd)
        if 'tempfiles' not in folders:
            sys.exit('No tempfiles/ found in directory, are you sure you`ve \
run pG-lt here?')
        avoid_pattern = "^\.|^log\.txt$|resolved_names|README|tempfiles"
        folders = [e for e in folders if not re.search(avoid_pattern, e)]
        folders = [os.path.join(wd, e) for e in folders]
        self.folders = folders
        self.backup_folders = folders
        self.parameter_keys = paradict.keys()
        self.parameter_keys.append('threads')
        self.gene_parameter_keys = genedict[genedict.keys()[0]].keys()
        self.gene_keys = genedict.keys()
        self.paradict = paradict
        self.genedict = genedict

    def _print(self, str):
        '''Special print'''
        if self.verbose:
            print(str)

    def _optionPrint(self, options):
        for option in options:
            self._print('    {0}'.format(option))

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

    def _resetstage(self, stage=None):
        """Reset `folders` to previous stage"""
        self._print('-'*70)
        while True:
            self._print('Reset `folders` to previous stages.')
            self._print('This will delete all files in folders from, but not \
including, the stage number given')
            if not stage:
                stage = raw_input('Enter stage (1-3): ')
            if int(stage) <= 3 and int(stage) > 0:
                break
            else:
                stage = None
                self._print('Invalid stage')
        counter = 0
        stagenames = {'1': '1_names', '2': '2_download', '3': '3_alignment',
                      '4': '4_phylogeny'}
        stages = range(int(stage)+1, 5)
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
        progress = self._readPickledFile(folder=self.wd, filename='progress.p')
        for s in stages:
            progress[str(s)] = 'not run'
        self._writePickledFile(folder=self.wd, filename='progress.p',
                               pickled=progress)
        self._print('    [{0}] folders reset to stage [{1}]'.format(counter,
                                                                    stage))

    def _resetparameters(self, key=None, value=None):
        '''Change key's value in all paradicts in `folders`'''
        self._print('-'*70)
        while True:
            self._print('Change parameters for all `folders`.\
Available options:')
            self._optionPrint(self.parameter_keys)
            if not key:
                key = raw_input('Enter parameter name of setting to change: ')
            if key not in self.parameter_keys:
                key = None
                self._print('Invalid parameter name!')
            else:
                break
        if not value:
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
        self._print('    [{0}] set to [{1}] for [{2}] folders'.
                    format(key, value, counter))

    def _resetgeneparameters(self, gene=None, key=None, value=None):
        '''Change key's value in a genedict in `folders`'''
        self._print('-'*70)
        while True:
            self._print('Enter name of gene to change. Available options:')
            self._optionPrint(self.gene_keys)
            if not gene:
                gene = raw_input('Enter name: ')
            if gene in self.gene_keys:
                break
            else:
                gene = None
                self._print('Invalid gene name')
        while True:
            self._print('Enter parameter name of setting to change. \
Available options:')
            self._optionPrint(self.gene_parameter_keys)
            if not key:
                key = raw_input('Enter name: ')
            if key in self.gene_parameter_keys:
                break
            else:
                key = None
                self._print('Invalid parameter name')
        if not value:
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
        self._print('    [{0}] set to [{1}] for [{2}] gene for [{3}] folders'.
                    format(key, value, gene, counter))

    def _setfolders(self, user_folders=None):
        '''Change `folders`'''
        self._print('-'*70)
        if not user_folders:
            user_folders = raw_input('Enter folder names separated by spaces: ')
        user_folders = user_folders.split(' ')
        backup_folders = [os.path.basename(e) for e in self.backup_folders]
        self.folders = []
        for folder in user_folders:
            if folder not in backup_folders:
                self._print('[{0}] not found -- did you spell it correctly?'.
                            format(folder))
            self.folders.append(os.path.join(self.wd, folder))

    def _restore(self, for_genedict=None, for_paradict=None, check=True):
        '''Save paradict and/or genedict in all `folders`'''
        while True:
            if for_genedict is None:
                for_genedict = raw_input('Restore default gene parameters? \
(y/n)')
            if for_genedict.lower() == 'y':
                for_genedict = True
                break
            elif for_genedict.lower() == 'n':
                for_genedict = False
                break
            else:
                self._print('Invalid argument')
        while True:
            if for_paradict is None:
                for_paradict = raw_input('Restore default parameters? \
(y/n)')
            if for_paradict.lower() == 'y':
                for_paradict = True
                break
            elif for_paradict.lower() == 'n':
                for_paradict = False
                break
            else:
                self._print('Invalid argument')
        if not for_genedict and not for_paradict:
            return(None)
        if check:
            _ = raw_input('Are you sure you want to restore to defaults? \
This is unreversable! Ctrl+c to exit, or hit enter to contiune')
            del _
        for folder in self.folders:
            if for_paradict:
                self._writePickledFile(folder=folder, filename='paradict.p',
                                       pickled=self.paradict)
            if for_genedict:
                self._writePickledFile(folder=folder, filename='genedict.p',
                                       pickled=self.genedict)

    def run(self):
        '''Run reset'''
        try:
            while True:
                self._print('-'*70)
                self._print('\n{0} RESET MODE {0}'.format(' '*29))
                self._print('-'*70)
                self._print(self.options_msg.format(self.wd))
                option = str(raw_input('Enter option number (1-5): '))
                if '1' == option:
                    self._resetstage()
                if '2' == option:
                    self._resetparameters()
                if '3' == option:
                    self._resetgeneparameters()
                if '4' == option:
                    self._setfolders()
                if '5' == option:
                    self._restore()
                # TODO: reset the ones that failed
        except KeyboardInterrupt:
            sys.exit('\nExiting reset mode ....')
