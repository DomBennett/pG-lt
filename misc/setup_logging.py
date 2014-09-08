import logging
logger = logging.getLogger('')
logger.setLevel(logging.INFO)
console = logging.StreamHandler()
console.setFormatter(logging.Formatter('%(message)s'))
logger.addHandler(console)

import mpe
mpe.stages.alignment_stage.run()