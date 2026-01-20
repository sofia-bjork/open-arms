#!/bin/bash -l

# Set the allocation to be charged for this job
# not required if you have set a default allocation
#SBATCH -A project-allocation

# The name of the script is myjob
#SBATCH -J openarms-miseq

# The partition
#SBATCH -p main

# 24 hours wall clock time will be given to this job
#SBATCH -t 16:00:00

# Number of nodes
#SBATCH --nodes=1

#SBATCH --mail-user=user@email.com
#SBATCH --mail-type=ALL

## Set the names for the error and output files. 
## It can be smart to set a path to these to your project directory, which you can do by adding that path right after the '=' sign
#SBATCH --error=/PATH/TO/open-arms/sbatch_jobs/job.%J.err
#SBATCH --output=/PATH/TO/open-arms/sbatch_jobs/job.%J.out

## module purge, then load necessary modules
## module purge 
## module load ...
## module load *apptainer module*


WORKDIR=/PATH/TO/open-arms;
TMPDIR=/PATH/TO/TEMP_DIR;
SANDBOX="$TMPDIR/temparms-sandbox"

apptainer exec --bind "$WORKDIR":/work --bind "$TMPDIR/conda_envs":/envs "$SANDBOX" bash /work/run/miseq_run.sh

