## MRes Project 2013
## Stage 3: Downloading sequences
## In: 5_phylodata | Out: 6_sumtrees
## 11/08/2013

## Print stage
print "\n\nThis is stage 3: download\n"

## Packages
import sys, os, re

## Dirs
input_dir = os.path.join(os.getcwd(), '6_screen')
output_dir = os.path.join(os.getcwd(), '7_sumtrees')
if not os.path.isdir(output_dir):
    os.mkdir(output_dir)
    
## Summarising
print "Looping through study phylogenies ..."
phylo_files = [p for p in os.listdir(input_dir) if re.search("\.tre", p)]
counter = 0
for phylo_file in phylo_files:
	out_phylo = re.sub("\.tre", "", phylo_file) + "_consensus"
	command = os.path.join(input_dir, phylo_file) + " > " + \
		os.path.join(output_dir, out_phylo)
	os.system("sumtrees.py --min-clade-freq=0.5 " + command)
	counter += 1
print "Done"
print "\n\nStage finished generated [{0}] consensus trees".format(counter)
