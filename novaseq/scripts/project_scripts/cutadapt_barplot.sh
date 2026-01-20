#!/usr/bin/env bash

OUTPUT_DIR="novaseq/COI/cutadapt_subset"
mkdir -p -v $OUTPUT_DIR

SUBSET_DIR=$OUTPUT_DIR/subset_gzfiles
mkdir -p -v $SUBSET_DIR

FASTQ_FILES="novaseq/COI/fastq_files"

FWD_SEQ=$FASTQ_FILES/DBQ_ABBPOSTA_2_1_HJ2LFDRX5.UDI172-BID05_clean.fastq.gz
REV_SEQ=$FASTQ_FILES/DBQ_ABBPOSTA_2_2_HJ2LFDRX5.UDI172-BID05_clean.fastq.gz

# zgrep -c ^@ $FWD_SEQ
# zgrep -c ^@ $REV_SEQ
# OUTPUT: 2 191 581
# 2 191 581 / 21 = 104 361
# zgrep -c $ $FWD_SEQ
# zgrep -c $ $REV_SEQ
# OUTPUT: 8 766 324
# 8 766 324 / 21 = 417 444
# 8 766 324 / 63 = 139 148 
# 63 * 4 = 252

# divide the 2 191 581 sequences (8 766 324 rows) into 63 files 
# (34 787 sequences and 139 148 rows in each file)
gzcat $FWD_SEQ | split -a 4 -l 139148 - $SUBSET_DIR/subset_ABBPOSTA_1.gz-
gzcat $REV_SEQ | split -a 4 -l 139148 - $SUBSET_DIR/subset_ABBPOSTA_2.gz-

# extract the subset 
suffixes=$(ls "$SUBSET_DIR"/subset_ABBPOSTA_1.gz-* | sed 's/.*-//')

# run cutadapt on each subset file
for suffix in $suffixes; do
    cutadapt --report=minimal -e 0.1 -g GGWACWGGWTGAACWGTWTAYCCYCC -G TANACYTCNGGRTGNCCRAARAAYCA \
             -o "$OUTPUT_DIR/${suffix}_ABBPOSTA_1.fastq" -p "$OUTPUT_DIR/${suffix}_ABBPOSTA_2.fastq" "$SUBSET_DIR/subset_ABBPOSTA_1.gz-${suffix}" "$SUBSET_DIR/subset_ABBPOSTA_2.gz-${suffix}" > "$OUTPUT_DIR/${suffix}_ABBPOSTA_cutadapt_summary.txt"
echo "${suffix} done"
done 

# remove .fastq files to save space
rm $OUTPUT_DIR/*.fastq
rm -r $SUBSET_DIR