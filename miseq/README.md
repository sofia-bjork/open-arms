# This is my project. I am doing a test run on the protocol from Nauras and then I will do the whole dataset. 

## These are the steps:

### Create conda environment from .yml file

`conda env create -f env/git_env.yml`

### Activate environment

`conda activate git_env`

### Download metadata files and specify parameters

`Rscript scripts/project_scripts/demultiplexed_summary.R -o Koster -p DMSO`

### Review and edit the extracted information, remove using grep or nano

`cat data/generated_meta/COI_demultiplexed_summary.csv`

### Create run number summary file 

`Rscript scripts/project_scripts/sequencing_run.R`

### Download the specified ENA fastq files

`Rscript scripts/original_scripts/ENADownload.R -f data/generated_meta/COI_ENA_accessions.txt -d COI/fastq_files`

`Rscript scripts/original_scripts/ENADownload_8digits.R -f data/generated_meta/COI_ENA_accessions_8digits.txt -d COI/fastq_files`

### Filter and trim using cutadapt 

`which swarm`

`Rscript scripts/project_scripts/COI_cutadapt_dada2.R -d CUTADAPT_DIR`

### Remove chimeras and singletons

`Rscript scripts/original_scripts/COI_chimera.R`

### Remove nuMTs using MACSE

`Rscript scripts/original_scripts/MACSE_align_pseudo.R -d PATH/TO/MACSE.jar </opt/miniconda3/envs/git_env/share/macse-2.07-0/macse_v2.07.jar>`

(macOS) `grep -w -A 1 -f COI/pseudo/nonpseudo.combined.names.txt COI/COI_nochim_nosingle_ASVs.fa > COI/pseudo/COI_nochim_nosingle_nopseudo.fa`

(Windows/GNU) `grep -w -A 1 -f COI/pseudo/nonpseudo.combined.names.txt COI/COI_nochim_nosingle_ASVs.fa > COI/pseudo/COI_nochim_nosingle_nopseudo.fa  --no-group-separator`

### Blank correction 

`Rscript scripts/original_scripts/COI_blank_corr.R`

### Clustering ASVs into MOTUs

`Rscript scripts/original_scripts/dereplication_headers.R`

`grep -v "^--" COI/blank_corr/COI_nochim_nosingle_nopseudo_nocontam.fa | awk 'NR%2==0' | paste -d'\n' COI/ASV_dereplicated.txt - > COI/COI_dereplicated_ASVs.fa`

`cd COI/`

`which swarm`

`PATH/TO/SWARM -d 13 -i swarm/internal.txt -o swarm/output.txt -s swarm/statistics.txt -u swarm/uclust.txt -w swarm/COI_cluster_reps.fa COI_dereplicated_ASVs.fa`

`cd ..`

### LULU curation

`Rscript scripts/original_scripts/MOTU_tables.R`

`awk -F'_' '{print $1}' COI/swarm/COI_cluster_reps.fa > COI/MOTU/COI_cluster_reps_lulu_ready.fa`

`makeblastdb -in COI/MOTU/COI_cluster_reps_lulu_ready.fa -parse_seqids -dbtype nucl`

`blastn -db COI/MOTU/COI_cluster_reps_lulu_ready.fa -outfmt '6 qseqid sseqid pident' -out COI/MOTU/match_list.txt -qcov_hsp_perc 80 -perc_identity 84 -query COI/MOTU/COI_cluster_reps_lulu_ready.fa`

`Rscript scripts/original_scripts/LULU_curation.R`

### Taxonomic assignment using BOLDigger3

`pip install boldigger3`

`pip install lxml_html_clean`

`scripts/project_scripts/BOLDigger_run.sh`

`Rscript scripts/project_scripts/BOLDigger_tax_table.R`
