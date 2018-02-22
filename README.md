Automagically process vector screening and file trimming locally.

Motivation
==========
* Do you find it annoying to run vecscreen and not be able to parse the results?
* Do you wish you could process the vector ID and then trim the sequences 5' or 3' from the vector/primer
* seqclean did not work to trim all of these out so I used vecscreen approach

This worked for me, happy to more generalize it if needed

1. Run BLAST jobs/vecscreen_blastn.sh - this script expects Trinity assemblies but can be changed
2. Run jobs/parse.sh - will parse the results and generated trimmed files for upload to TSA or other GenBank db

AUTHOR
======
Jason Stajich - jason.stajich[AT]ucr.edu
