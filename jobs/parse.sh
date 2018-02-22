#!/bin/bash
#SBATCH --ntasks 6 --nodes 1 --mem 8G --out parse.log -J cleanvectparse  -p short 

# fix this if installed do not need a relative folder
parallel -j 6 ./scripts/clean_vect_blastn.py {} {}.vecscreen.blastn ::: *Trinity_Assembly.fasta
