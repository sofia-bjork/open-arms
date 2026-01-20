#!/bin/bash

### specify path to BOLDigger output directory and create if not present ###
BOLDIGGER_DIR="COI/BOLDigger"
mkdir -p -v $BOLDIGGER_DIR

### specify path to LULU curation output directory ###
LULU_DIR="COI/MOTU"

### generate FASTA file with MOTUs remaining after LULU curation ###
grep -w -A 1 -f $LULU_DIR/lulu_curated_headers.txt $LULU_DIR/COI_cluster_reps_lulu_ready.fa | \
grep -v "^--" > $BOLDIGGER_DIR/COI_cluster_reps_lulu_curated.fa 

### run BOLDigger on exhaustive search (--mode 3) and match to public database (--db 1) ###
boldigger3 identify $BOLDIGGER_DIR/COI_cluster_reps_lulu_curated.fa --db 1 --mode 3

#################################################################################################

### specify path to BOLDigger output directory and create if not present ###
BOLDIGGER_PRIV="COI/BOLDigger_private"
mkdir -p -v $BOLDIGGER_PRIV

cp $BOLDIGGER_DIR/COI_cluster_reps_lulu_curated.fa  $BOLDIGGER_PRIV/COI_cluster_reps_lulu_curated.fa

### run BOLDigger on exhaustive search (--mode 3) and match to private database (--db 2) ###
boldigger3 identify $BOLDIGGER_DIR/COI_cluster_reps_lulu_curated.fa --db 2 --mode 3
