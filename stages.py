import binaries

STAGES = { '1'   : (binaries.python27, '-u', '1_taxids.py'),
           '2'   : (binaries.python27, '-u', '2_download.py'),
           '3'   : (binaries.python27, '-u', '3_alignments.py'),
           '4'   : (binaries.python27, '-u', '4_phylogenies.py')
         }