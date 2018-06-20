#!/usr/bin/env python3
import sys, csv, re

from Bio import SeqIO
from Bio.Seq import Seq

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

                cols = exline.split(maxsplit=4)

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
            print(record.id, "before is ",len(record), "long")
            if len(trimloc) > 1:
                print("more than one trim location ... maybe a chimera?",record.id)
            else:
                locs = trimloc[0].split("..")
                left = int( locs[0] ) - 1
                right = int( locs[0] )
                if left == 0:
                    record = record[right-1:]
                elif right == len(record):
                    record = record[:left]
                else:
                    # internal slicing
                    temprecord = Seq(record[:left] + record[right-1:],DNAAlphabet())
                    temprecord.id = record.id
                    print(record.id, len(record))
                    record = temprecord
                    print(len(record))
            print(record.id, "after is ",len(record), "long")
        elif record.id in duplicates:
           log.write("Skipping %s as is considered a duplicate\n" % record.id)
           continue

        SeqIO.write(record, output_handle, "fasta")
