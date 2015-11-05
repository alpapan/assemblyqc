#!/bin/bash
bunzip2 -k ../db/archae_bacteria_plasmid_viral/taxonomy/*bz2 2>/dev/null

bunzip2 -k ../Ccap01172013-genome.fa.bz2 2>/dev/null

check_contamination.pl -i Ccap01172013-genome.fa -d ../kraken/archae_bacteria_plasmid_viral
