#!/usr/bin/env python3
import sys, csv, re, operator
# biopython needed
from Bio import SeqIO

# this code is based on this post on stackexchange
# https://codereview.stackexchange.com/questions/69242/merging-overlapping-intervals
def merge_intervals(intervals):
    """
    A simple algorithm can be used:
    1. Sort the intervals in increasing order
    2. Push the first interval on the stack
    3. Iterate through intervals and for each one compare current interval
       with the top of the stack and:
       A. If current interval does not overlap, push on to stack
       B. If current interval does overlap, merge both intervals in to one
          and push on to stack
    4. At the end return stack
    """
    merged = []
    sorted_by_lower_bound = sorted(intervals, key=operator.itemgetter(0))

    if not sorted_by_lower_bound:  # no intervals to merge
        return

    for higher in sorted_by_lower_bound:
        if not merged:
            merged.append(higher)
        else:
            lower = merged[-1]
            #added this condition branch
            if higher[0] - lower[1] == 1:
                merged[-1] = (lower[0], higher[1])  # replace by merged interval
            #end addition. Also changed if below to elif
            # test for intersection between lower and higher:
            # we know via sorting that lower[0] <= higher[0]
            elif higher[0] <= lower[1]:
                upper_bound = max(lower[1], higher[1])
                merged[-1] = (lower[0], upper_bound)  # replace by merged interval
            else:
                merged.append(higher)
    return merged

fastafile = sys.argv[1]
blastn = sys.argv[2]
cleaned = fastafile + ".clean.fsa"
logging = fastafile + ".parse.log"

excludes = {}

with open(blastn,"r") as vectab:
    rdr = csv.reader(vectab,delimiter="\t")
    for row in rdr:
        idin = row[0]
        loc = [int(row[8]), int(row[9])]
        if loc[0] > loc[1]:
            loc = [loc[1],loc[0]]

        if idin not in excludes:
            excludes[idin] = []

        excludes[idin].append(loc)

with open(cleaned, "w") as output_handle, open(logging,"w") as log:
    for record in SeqIO.parse(fastafile, "fasta"):
        record.description = ""
        if record.id in excludes:
            trimloc = excludes[record.id]

            if len(trimloc) > 1:
                trimloc = merge_intervals(trimloc)

            seqlen = len(record)
            for l in trimloc:
                if l[0] == 1:
                    # 5' trimming
                    record = record[l[1]:]
                elif l[1] == seqlen :
                    # 3' trimming
                    record = record[:l[0]]
                else:
                    if( l[0]- 1  < seqlen - l[1]):
                        record = record[l[1]:]
                        log.write( "%s [%d,%d] trimming 5' though does not match 1 or len=%d\n" % 
                                   (record.id,l[0],l[1],seqlen))
                    else:
                        log.write( "%s [%d,%d] trimming 3' though does not match 1 or len=%d\n" % 
                                   (record.id,l[0],l[1],seqlen))
                        record = record[:l[0]]

        if(len(record) >= 200):
            SeqIO.write(record, output_handle, "fasta")
