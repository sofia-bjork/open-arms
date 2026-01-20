#!/usr/bin/env bash

export PATH=/envs/git_env/bin:$PATH

set -e

cd /PATH/TO/open-arms

# download metadata files and specify parameters of interest
# -o observatory -p preservative -u UnitID -f fraction
Rscript miseq/scripts/project_scripts/demultiplexed_summary.R -o Koster -p DMSO

# create run number summary file
Rscript scripts/project_scripts/sequencing_run.R

# download the specified ENA fastq files
Rscript miseq/scripts/original_scripts/ENADownload.R -f metadata/generated_meta/COI_ENA_accessions.txt -d miseq/COI/fastq_files
Rscript miseq/scripts/original_scripts/ENADownload_8digits.R -f metadata/generated_meta/COI_ENA_accessions_8digits.txt -d miseq/COI/fastq_files

# filter and trim using cutadapt and dada2
Rscript miseq/scripts/project_scripts/COI_cutadapt_dada2.R -d /opt/miniconda3/envs/git_env/bin/cutadapt

# remove chimeras and singletons
Rscript miseq/scripts/original_scripts/COI_chimera.R

# identify and remove nuclear mitochondrial DNA pseudogenes (nuMTs)
Rscript miseq/scripts/original_scripts/MACSE_align_pseudo.R -d /opt/miniconda3/envs/git_env/share/macse-2.07-0/macse_v2.07.jar

# subset non-nuMT ASVs from the COI_nochim_nosingle.fa file
grep -w -A 1 -f miseq/COI/pseudo/nonpseudo.combined.names.txt miseq/COI/COI_nochim_nosingle_ASVs.fa > miseq/COI/pseudo/COI_nochim_nosingle_nopseudo.fa  --no-group-separator

# negative control correction
Rscript miseq/scripts/original_scripts/COI_blank_corr.R

# subset the ASVs remaining after negative control correction
grep -w -A 1 -f miseq/COI/blank_corr/no_contam_headers_COI.txt miseq/COI/pseudo/COI_nochim_nosingle_nopseudo.fa  --no-group-separator > miseq/COI/blank_corr/COI_nochim_nosingle_nopseudo_nocontam.fa 

# generate headers with the abundance of each ASV included
Rscript miseq/scripts/original_scripts/dereplication_headers.R

# replace headers with dereplicated header names
grep -v "^--" miseq/COI/blank_corr/COI_nochim_nosingle_nopseudo_nocontam.fa | awk 'NR%2==0' | paste -d'\n' miseq/COI/ASV_dereplicated.txt - > miseq/COI/COI_dereplicated_ASVs.fa

# change directory for swarm to run 
pushd miseq/COI/

# cluster ASVs into MOTUs using swarm
/opt/miniconda3/envs/main_env/bin/swarm -d 13 -i swarm/internal.txt -o swarm/output.txt -s swarm/statistics.txt -u swarm/uclust.txt -w swarm/COI_cluster_reps.fa COI_dereplicated_ASVs.fa

# move back to root
popd

# generate MOTU tables from swarm output 
Rscript miseq/scripts/original_scripts/MOTU_tables.R

# remove read abundance line from sequence header
awk -F'_' '{print $1}' miseq/COI/swarm/COI_cluster_reps.fa > miseq/COI/MOTU/COI_cluster_reps_lulu_ready.fa

# generate match lists using BLASTn
makeblastdb -in miseq/COI/MOTU/COI_cluster_reps_lulu_ready.fa -parse_seqids -dbtype nucl
blastn -db miseq/COI/MOTU/COI_cluster_reps_lulu_ready.fa -outfmt '6 qseqid sseqid pident' -out miseq/COI/MOTU/match_list.txt -qcov_hsp_perc 80 -perc_identity 84 -query miseq/COI/MOTU/COI_cluster_reps_lulu_ready.fa

# LULU curation
Rscript miseq/scripts/original_scripts/LULU_curation.R

# generate fasta files with remaining MOTUs after LULU curation
grep -w -A 1 -f miseq/COI/MOTU/lulu_curated_headers.txt miseq/COI/MOTU/COI_cluster_reps_lulu_ready.fa --no-group-separator > miseq/COI/MOTU/COI_cluster_reps_lulu_curated.fa

# taxonomic assignment using BOLDigger3 on public animal library (--db 1) on exhaustive search mode (--mode 3)
pip install boldigger3==2.1.4
pip install lxml-html-clean==0.4.3
boldigger3 identify miseq/COI/MOTU/COI_cluster_reps_lulu_curated.fa --db 1 --mode 3

# create taxonomic table from boldigger output 
Rscript miseq/scripts/project_scripts/BOLDigger_tax_table.R