#!/bin/bash
#SBATCH --nodes 1 --ntasks 16 --mem 8G --out vecscreen.log -p short --time 2:00:00

module load ncbi-blast

# NCBI parameters for vector screen, doing this with parallel for now

JOBS=4
CPUS_PER_JOB=4
DB=UniVec

parallel -j $JOBS blastn -task blastn -reward 1 -penalty -5 -gapopen 3 -gapextend 3 -dust yes -soft_masking true -evalue 700 -searchsp 1750000000000 \
 -db $DB -outfmt 6 -num_threads $CPUS_PER_JOB -query {} -out {}.vecscreen.blastn ::: *Trinity_Assembly.fasta
