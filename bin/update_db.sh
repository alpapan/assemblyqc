#!/bin/bash

# NB Requires jellyfish 1
# aptitude install jellyfish

NCPUS=15
FILE=new_db_for_kraken

# these are the commands, change the relevant fasta file
/usr/bin/jellyfish count -m 31  -C -t $NCPUS -s 4G -o $FILE.jdb $FILE.fsa
/usr/bin/jellyfish merge --output=$FILE.jdb $FILE.jdb_*

../kraken/db_sort -z -t $NCPUS -n 15 -d $FILE.jdb \
 -o ../db/archae_bacteria_plasmid_viral/taxonomy/database.kdb -i ../db/archae_bacteria_plasmid_viral/taxonomy/database.idx &

../kraken/report_gi_numbers.pl $FILE.fsa > $FILE.gis
../kraken/make_seqid_to_taxid_map ../db/archae_bacteria_plasmid_viral/taxonomy/gi_taxid_nucl.dmp $FILE.gis > $FILE.map


../kraken/set_lcas -x -d ../db/archae_bacteria_plasmid_viral/database.kdb -i ../db/archae_bacteria_plasmid_viral/taxonomy/database.idx \
 -n ../db/archae_bacteria_plasmid_viral/taxonomy/nodes.dmp -t $NCPUS -m $FILE.map -F $FILE.fsa
