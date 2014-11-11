# Set up runner

import pglt
import os


Runner = pglt.tools.system_tools.Runner
folders = ['SE1_2007__Parra_1']
runner = Runner(folders, 1, os.getcwd())
runner.run(folders, ['1', '2'])
