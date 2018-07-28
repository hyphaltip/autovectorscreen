#!/usr/bin/env python3
import sys, csv, re, operator

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

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
ncbirpt = sys.argv[2]
cleaned = fastafile + ".filtered_clean.fsa"
logging = fastafile + ".parse_tagged.log"

excludes = {}
trims = {}
duplicates = {}

with open(ncbirpt,"r") as rpt:    
    status = ""
    for line in rpt:
        if  line.startswith("Exclude:"):
            status = "E"
            header = rpt.readline()
            for exline in rpt:
                if exline.startswith("\n"):
                    break # end the exclude loop

                cols = exline.split(maxsplit=3)
                excludes[cols[0]] = 1

        elif line.startswith("Trim:"):
            status = "T"
            header = rpt.readline()
            for tline in rpt:
                if tline.startswith("\n"):
                    break # end the trim

                cols = tline.split(maxsplit=4)
                if cols[0] not in trims:
                    trims[cols[0]] = []
                locs = cols[2].split(",")
                #print("locs are %s" % locs)
                for l in locs:
                    loc1 = [ int(x) for x in l.split("..")]
                    #print ("adding %s to list for location %s" % ( cols[0],loc1))
                    trims[cols[0]].append(loc1)

        elif line.startswith("Duplicated:"):
            status = "D"
            header = rpt.readline()
            for dline in rpt:
                if dline.startswith("\n"):
                    break # end the trim
                
                cols = dline.split()
                len = cols.pop()
                len = cols.pop()
                for col in cols[1:]:
                    col = re.sub("lcl\|","",col)
#                    print("line is ",dline," marking ",col,
#                          " to remove as dup")
                    duplicates[col] = 1
                
with open(cleaned, "w") as output_handle, open(logging,"w") as log:
    for record in SeqIO.parse(fastafile, "fasta"):
        record.description = ""
        if record.id in excludes:
            #skip this record
            log.write("Skipping %s as is considered a contaminant\n" % record.id)
            continue
        elif record.id in trims:
            trimloc = trims[record.id]
            print(record.id, "Begin trim loop for %s is len %d" %(record.id,len(record)))
            if len(trimloc) > 1:
                trimloc = sorted(merge_intervals(trimloc),
                                 reverse=True,
                                 key=lambda locitem: locitem[0])
#                print("trim loc is ",trimloc)
#                print("more than one trim location ... maybe a chimera?",record.id)
            for loc in trimloc:
                left = int( loc[0] ) - 1
                right = int( loc[1] )
                newrecord = Seq("",generic_dna)
                print("trimming %d to %d in %s len=%d" 
                      % (left,right,record.id,len(record)))
                if left == 0:
                    newrecord = record[right-1:]
                elif right == len(record):
                    newrecord = record[:left]
                else:
                    # internal slicing
#                    print("-->internal slicing :%d .. %d:" % (left,right-1))
#                    print('left string is', record[0:left])
#                    print('right string is', record[right-1:])
                    newrecord = record[0:left] + record[right-1:]
                    newrecord.id = record.id
#                    print(" -- new len is",newrecord.id, len(newrecord),
#                          newrecord)
                record = newrecord
#                print(record.id, "after is ",len(record), "long")
        
            print("done with trimming loop length %s is now %d"
                  % (record.id, len(record)) )
                                                                   
        elif record.id in duplicates:
            log.write("Skipping %s as is considered a duplicate\n" % record.id)
            continue

        if(len(record) >= 200):
            SeqIO.write(record, output_handle, "fasta")

