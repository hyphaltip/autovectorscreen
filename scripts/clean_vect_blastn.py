#!/usr/bin/env python3
import sys, csv, re, operator,os
# biopython needed
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

percent_id_cutoff = 98
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
prefix=os.path.splitext(fastafile)[0]
cleaned = prefix + ".clean.fsa"
logging = prefix + ".parse.log"

excludes = {}

with open(blastn,"r") as vectab:
    rdr = csv.reader(vectab,delimiter="\t")
    for row in rdr:
        idin = row[0]
        pid = row[2]
        if float(pid) < percent_id_cutoff:
            continue
            # skip this line since PID is lower than needed

        loc = [int(row[6]), int(row[7])]
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
                trimloc = sorted(merge_intervals(trimloc),
                                 reverse=True,
                                 key=lambda locitem: locitem[0])

            seqlen = len(record)
            for loc in trimloc:
                left = int( loc[0] ) - 1
                right = int( loc[1] )
                newrecord = Seq("",generic_dna)
                log.write("trimming %d to %d in %s len=%d" 
                      % (left,right,record.id,len(record)))
                if left == 0:
                    newrecord = record[right-1:]
                elif right == len(record):
                    newrecord = record[:left]
                else:
                    # internal slicing
                    log.write("-->internal slicing :%d .. %d:" % (left,right-1))
                    log.write('  left string is %s' % record[0:left])
                    log.write('  right string is %s' % record[right-1:])
                    newrecord = record[0:left] + record[right-1:]
                    newrecord.id = record.id
                    log.write(" -- new len for %s is %d: %s" % (newrecord.id, len(newrecord),newrecord))
                record = newrecord

        if(len(record) >= 200):
            SeqIO.write(record, output_handle, "fasta")
