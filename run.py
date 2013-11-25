""" Usage: run.py stage

"""

import os
import platform
import socket
import stat
import subprocess
import sys

import stages

from datetime import datetime

class StageError(Exception):
    pass

class Tee(object):
    # http://stackoverflow.com/questions/616645/how-do-i-duplicate-sys-stdout-to-a-log-file-in-python/616686#616686
    def __init__(self, path):
        self.file = open(path, 'w')
        self.stdout,sys.stdout = sys.stdout,self
        self.stderr,sys.stderr = sys.stderr,self

    def __enter__(self):
        pass

    def __exit__(self, type, value, traceback):
        sys.stdout,sys.stderr = self.stdout,self.stderr
        self.file.close()

    def write(self, data):
        self.file.write(data)
        self.file.flush()
        self.stdout.write(data)
        self.stdout.flush()

class Stage(object):
    STAGES = stages.STAGES

    def __init__(self, stage):
        if stage not in self.STAGES:
            raise StageError('Stage [{0}] not recognised'.format(stage))
        else:
            self.stage = stage
            self.output_dir = os.path.splitext(self.STAGES[stage][-1])[0]

    def _prime(self):
        # Creates self.output_dir. An error is raised if it already exists.
        if os.path.isdir(self.output_dir):
            raise StageError('Output directory [' + self.output_dir + '] ' + 
                             'already exists. Has this stage already been run?')
        else:
            os.mkdir(self.output_dir)
            return Tee(os.path.join(self.output_dir, 'log.txt'))

    def _time_string(self):
        return datetime.today().strftime("%A, %d %B %Y %I:%M%p")

    def _finished(self):
        # Make output files read only
        for root,dirs,files in os.walk(self.output_dir):
            for f in files:
                os.chmod(os.path.join(root, f), stat.S_IREAD)

    def _cmd(self, args):
        s = subprocess.Popen(args, stdout=subprocess.PIPE, 
                             stderr=subprocess.STDOUT)
        while True:
            line = s.stdout.readline()
            exitcode = s.poll()
            line = line[:-1]
            if (not line) and (exitcode is not None):
                break
            elif line:
                print line

        if exitcode:
            sys.exit(exitcode)

    def run(self):
        with self._prime():
            print 'Stage [{0}] started at [{1}]'.format(self.stage, 
                                                        self._time_string())
            print 'Running on [{0}] [{1}]'.format(platform.node(), 
                                                  platform.platform())
            args = self.STAGES[self.stage]
            self._cmd( [args[0], '--version'] )
            self._cmd(args)
            print 'Stage finished at [{0}]'.format(self._time_string())

        self._finished()

    @classmethod
    def run_all(klass):
        print 'Running all stages'
        for s in sorted(Stage.STAGES.keys()):
            Stage(s).run()

def main():
    if 2!=len(sys.argv):
        print __doc__
        sys.exit(1)
    elif 'all'==sys.argv[1].lower():
        Stage.run_all()
    else:
        Stage(sys.argv[1]).run()

if __name__=='__main__':
    main()

