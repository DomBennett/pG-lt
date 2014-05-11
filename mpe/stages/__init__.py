import os
from mpe import _ROOT

path = os.path.join(_ROOT, 'stages')
STAGES = {'1':os.path.join(path, '1_names.py'),
'2':os.path.join(path, '2_download.py'),
'3':os.path.join(path,'3_alignment.py'),
'4':os.path.join(path,'phylogeny.py')}

del os
del path
del _ROOT