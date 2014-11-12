# Set up runner

import os
from pglt.tools.setup_tools import getDirs
from pglt.tools.setup_tools import setUpLogging
from pglt.tools.system_tools import Runner

base_logger = setUpLogging(False, False)
folders = getDirs(base_logger)
runner = Runner(folders, 1, os.getcwd(), 'dominic.john.bennett@gmail.com')
runner.setup(folders)
runner.run(folders, ['1'])


threads = []
for i in range(nworkers):
    t = threading.Thread(target=runner._worker)
    threads.append(t)
    t.daemon = True
    t.start()
