#!/bin/bash
#SBATCH --nodes 1 --ntasks 16 --mem 8G --out vecscreen.log -p short --time 2:00:00

module load parallel
module load ncbi-blast
module unload python
module unload python
module load python/3

# NCBI parameters for vector screen, doing this with parallel for now

JOBS=2
CPUS_PER_JOB=8
DB=/srv/projects/db/UniVec/UniVec
echo $0
SCRIPTS=$(dirname `dirname $0`)/scripts
echo $SCRIPTS

PID_CUTOFF=98
blast_n_parse() {
    PREF=$(basename $1 .fasta)
    ROUND=0
    FASTA=$PREF.r$ROUND.fna
    ln -s $1 $FASTA

    RUN=1
    while [[ $RUN == 1 ]];
    do	
	echo "RUN=$RUN ROUND=$ROUND PREF=$PREF FASTA=$FASTA"
	if [ ! -f $PREF.r$ROUND.vecscreen.tab ]; then
	    blastn -task blastn -reward 1 -penalty -5 -gapopen 3 \
		-gapextend 3 -dust yes -soft_masking true -evalue 700 -searchsp 1750000000000 \
		-db $DB -outfmt 6 -num_threads $CPUS_PER_JOB -query $FASTA \
		-out $PREF.r$ROUND.vecscreen.tab
	fi

	COUNT=0
	if [ -f $PREF.r$ROUND.vecscreen.tab ]; then
	    COUNT=$(awk -v cutoff=$PID_CUTOFF 'BEGIN { sum = 0 } { if ( $3 >= cutoff ) sum +=1 } END { print sum }' $PREF.r$ROUND.vecscreen.tab)
	fi

	if [[ $COUNT == 0 ]]; then
	    RUN=0
	else
	    NEXTROUND=$(expr $ROUND + 1)
	    $SCRIPTS/clean_vect_blastn.py $FASTA $PREF.r$ROUND.vecscreen.tab
	    ln -s $PREF.r$ROUND.clean.fsa $PREF.r$NEXTROUND.fna
	    ROUND=$NEXTROUND
	    FASTA=$PREF.r$ROUND.fna
	fi
    done
}
export -f blast_n_parse

echo "ready to run"
#parallel -j 1 blast_n_parse ::: *.fasta
for file in *.fasta
do
    echo $file
    blast_n_parse $file
done

# fix this if installed do not need a relative folder
#parallel -j 6 ./scripts/clean_vect_blastn.py {} {}.vecscreen.blastn ::: *.fasta
