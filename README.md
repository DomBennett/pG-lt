# pG-lt
[![Coverage Status](https://coveralls.io/repos/DomBennett/pG-lt/badge.svg?branch=master)](https://coveralls.io/r/DomBennett/pG-lt?branch=master)

'[phyloGenerator][pG]-light' or `pG-lt` is a python package for the automated
generation of phylogenies from multiple lists of taxon names.

The package implements a pipeline that uses 'mass phylogeny estimation' to
generate distributions of phylogenies from multiple sequences and alignments
with minimum human-input.

The pipeline consists of four stages: names resolution, sequence download,
sequence alignment and phylogeny estimation.

pG-lt is a member of the phyloGenerator family.

## Installation and usage

Click [here][wiki], to see the wiki for information on how to use and install.

## Features

* One-line implementation: `run_pglt.py -e [EMAIL ADDRESS]`
* Automatic or user specified parameters
* Fuzzy taxonomic name matching
* Automatic outgroup selection
* Heuristic gene selection
* Automatic sequence download with taxonomic context
* Heuristic alignment and phylogeny checking
* Rate partitioning between genes and codon positions
* HPC compatible
* Runs on Windows, OSX and Linux

## External dependencies

* Alignment tools: [MAFFT][mafft] and stand-alone [BLAST][blast]
* Phylogeny tools: [RAxML][raxml]

## Version

*pG-lt is still in its developmental stage.*

## License

GPL v.2

## Authors
[Dom Bennett][db], [Will Pearse][wp] and [Lawrence Hudson][lh]

<!-- References -->
[db]: https://github.com/DomBennett
[wp]: https://github.com/willpearse
[lh]: https://github.com/quicklizard99
[pG]: http://willpearse.github.io/phyloGenerator/
[wiki]: https://github.com/DomBennett/pG-lt/wiki
[mafft]: http://mafft.cbrc.jp/alignment/software/
[raxml]: https://github.com/stamatak/standard-RAxML
[blast]: http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download
