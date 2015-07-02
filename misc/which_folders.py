# 03/07/2015
# Dom Bennett
# Which folders succeeded or failed for which stages

# LIBS
import os
import sys
import sys

# FUNCTIONS
def checkProgress(folder, stage, succeed):
    """Return T/F for progress of folder given stage and succeed"""
    filepath = os.path.join(folder, 'tempfiles', 'progress.p')
    if not os.path.isfile(filepath):
        return False
    with open(filepath, "rb") as file:
        progress = pickle.load(file)
    if progress[stage] == 'not run':
        return False
    if progress[stage] == 'success':
        res = True
    else:
        res = False
    if not succeed:
        res = not res
    return res

# MAIN
if __name__ == '__main__':
    stage = sys.argv[1]
    succeed = sys.argv[2]
    if succeed.lower() == 'y':
        succeed = True
    elif succeed.lower() == 'n':
        succeed = False
    else:
        sys.exit('Second argumnet must be: y or n')
    folders = os.listdir('.')
    counter = 0
    for folder in folders:
        if checkProgress(folder, stage, succeed):
            print folder
            counter += 1
    sys.exit('Complete. [' + str (counter), '] folders.')
