#! /bin/usr/env python
# D.J. Bennett
# 24/03/2014
"""
pglt system tools
"""

# PACKAGES
import subprocess
import threading
import os
import Queue
from datetime import datetime
from setup_tools import setUpLogging
from setup_tools import tearDownLogging

# MESSAGES
priming_msg = '\nERROR: The program was unable to start due to a \
problem with the files and folders in the study directory. Check the \
parameters and gene parameters .csv for any potential conflicts.'
toofewspecies_msg = '\nERROR: The program halted as there are too few \
species left of phylogeny building -- five is the minimum. You may \
have started with too few names, or names given could not be \
taxonomically resolved or there may be too little sequence data \
available.'
taxonomicrank_msg = '\nERROR: It is likely that one or more names\
have been resolved incorrectly, as such the parent taxonomic group \
has been set to Eukaryotes which is too high a taxonomic rank for \
phylogenetic analysis. Consider adding a parent ID to the \
parameters.csv to prevent incorrect names resolution or reducing the \
taxonomic diversity of the analysis names.'
outgroup_msg = '\nERROR: The outgroup has been dropped. This may be \
due to too few sequence data available for outgroup or a failure to \
align sequences that are available. If outgroup has been \
automatically selected, consider manually choosing an outgroup.'
raxml_msg = '\nERROR: Generated maxtrys poor phylogenies \
consecutively, consider reducing maxrttsd.'
unexpected_msg = '\nERROR: The following unexpected error occurred:\n\
\"{0}\" \n\
Please email details to the program maintainer for help.'


# ERROR CLASSES
class StageError(Exception):
    pass


class TooFewSpeciesError(Exception):
    pass


class TaxonomicRankError(Exception):
    pass


class PrimingError(Exception):
    pass


class OutgroupError(Exception):
    pass


class MafftError(Exception):
    pass


class RAxMLError(Exception):
    pass


class TrysError(Exception):
    pass


# OTHER CLASSES
class Stager(object):
    """Stager class : runs each file in stage folder. Adapted from\
 code written by L. Hudson."""

    def __init__(self, wd, stage, verbose=False, debug=False):
        # STAGES is added to Stager at __init__.py
        if stage not in self.STAGES:
            raise StageError('Stage [{0}] not recognised'.format(stage))
        else:
            self.wd = wd
            self.folder = os.path.split(wd)[-1]
            self.stage = stage
            # dir is second element of tuple
            self.output_dir = os.path.join(wd, self.STAGES[stage][1])
            self.verbose = verbose
            self.debug = debug

    def _start(self):
        self.logger.info('-' * 70)
        self.logger.info('Stage [{0}] started at [{1}]'.
                         format(self.stage, self._time_string()))
        self.logger.info('-' * 70)

    def _end(self):
        self.logger.info('-' * 70)
        self.logger.info('Stage [{0}] finished at [{1}]'.
                         format(self.stage, self._time_string()))
        self.logger.info('-' * 70 + '\n\n')

    def _error(self, msg):
        """Return true when error raised, log informative message"""
        self.logger.error(msg)
        self.logger.info('.... Moving to next folder')
        self.logger.info('Stage [{0}] unfinished at [{1}]'.
                         format(self.stage, self._time_string()))
        self.logger.info('-' * 70 + '\n\n')
        return True

    def _time_string(self):
        return datetime.today().strftime("%A, %d %B %Y %I:%M%p")

    def _cmd(self):
        """Run stage command. Catch errors raised."""
        error_raised = None
        try:
            # function is first element of tuple
            # pass wd and logger
            self.STAGES[self.stage][0](self.wd, self.logger)
        except PrimingError:
            error_raised = self._error(priming_msg)
        except TooFewSpeciesError:
            error_raised = self._error(toofewspecies_msg)
        except TaxonomicRankError:
            error_raised = self._error(taxonomicrank_msg)
        except OutgroupError:
            error_raised = self._error(outgroup_msg)
        except RAxMLError:
            error_raised = self._error(raxml_msg)
        except Exception as unexpected_error:
            error_raised = self._error(unexpected_msg.format(unexpected_error))
        return error_raised

    def run(self):
        # make sure dir exists
        if not os.path.isdir(self.output_dir):
            os.mkdir(self.output_dir)
        # set up a logger
        self.logger = setUpLogging(verbose=self.verbose, debug=self.debug,
                                   logname=self.folder,
                                   directory=self.output_dir)
        # log system info
        self._start()
        # run stage
        failed = self._cmd()
        if not failed:
            # log end time
            self._end()
        # remove logger
        tearDownLogging(self.folder)
        return failed

    @classmethod
    def run_all(klass, wd, stage, verbose):
        for s in sorted(Stager.STAGES.keys()[stage:]):
            if not verbose:
                print '  Stage [{0}]'.format(Stager.STAGES[s][1])
            Stager(wd, s).run()


class Runner(object):
    """Runner class : run stages across folders"""
    def __init__(self, folders, nworkers, wd):
        self.wd = wd
        self.nworkers = nworkers
        self.folders = folders
        self.q = Queue.Queue(maxsize=nworkers)

    def _worker(self):
        while True:
            # get folder and stages from queue
            folder, stages = self.q.get()
            # get a working dir for folder
            stage_wd = os.path.join(self.wd, folder)
            # run each stage for folder
            for stage in stages:
                stager = Stager(stage_wd, stage)
                failed = stager.run()
                # if failed, remove folder from list
                if failed:
                    self.folders.remove(folder)
                    break
            self.q.task_done()

    def setup(self, folders):
        """Setup files across folders"""
        pass

    def run(self, folders, stages, parallel=False):
        """Run folders and stages"""
        if parallel:
            nworkers = self.nworkers
        else:
            nworkers = 1
        # create nworkers workers
        for i in range(nworkers):
            t = threading.Thread(target=self._worker)
            t.daemon = True
            t.start()
        # set workers running across all folders for stages
        for folder in folders:
            self.q.put((folder, stages))


class TerminationPipe(object):
    """TerminationPipe class : exectute background programs. Adapted pG code \
written by W.D. Pearse."""
    def __init__(self, cmd, timeout=99999, silent=True):
        self.cmd = cmd
        self.timeout = timeout
        self.process = None
        self.output = None
        self.failure = False
        self.stderr = 'EMPTY'
        self.stdout = 'EMPTY'
        self.silent = silent

    def run(self):
        def silentTarget():
            self.process = subprocess.Popen(self.cmd, stdout=subprocess.PIPE,
                                            shell=True, stderr=subprocess.PIPE)
            self.output = self.process.communicate()

        def loudTarget():
            self.process = subprocess.Popen(self.cmd, shell=False)
            self.output = self.process.communicate()
        if self.silent:
            thread = threading.Thread(target=silentTarget)
        else:
            thread = threading.Thread(target=loudTarget)
        thread.start()
        # TODO: how to handle the main thread being killed?
        thread.join(self.timeout)
        if thread.is_alive():
            self.process.terminate()
            thread.join()
            self.failure = True