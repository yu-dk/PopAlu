#!/usr/bin/python

import random
from sys import argv

# -----------------------------------------------------------------------------
# initialize parameters
# -----------------------------------------------------------------------------

if len(argv) != 3 and len(argv) != 4:
   print "Usage:", argv[0], "<ins.fa> <hap.vcf> [<seed>]"
   exit()

# -----------------------------------------------------------------------------

insfile = argv[1]
hapfile = argv[2]

seed = 0
if len(argv) > 3:
    seed = int(argv[3])

random.seed(seed)

# -----------------------------------------------------------------------------
# iterate insfile and append vcf records to hapfile
# -----------------------------------------------------------------------------

with open(hapfile, "a") as outfile:
    with open(insfile, "r") as infile:
        for line in infile:
            # skip comment lines starting with '#'
            if line[0] == '#':
                continue
            
            # split line into fields
            chr, pos, id, ref, alt, qual, filter, info, format, sim = line.split()
            
            # split the INFO field and extract allele frequency
            infos = info.split(";")
            i = 0
            while i < len(infos) and infos[i][:3] != "AF=":
                i += 1
            if i == len(infos):
                print "ERROR: Could not find allele frequency in INFO field."
                exit(1)
            freq = float(infos[i][3:])

            # append insertion to vcf with probability proportional to frequency
            if random.random() <= freq:
                outfile.write(chr + "\t" + pos + "\t" + id + "\t" + ref + "\t" + alt + "\t")
                outfile.write(qual + "\t" + filter + "\t" + info + "\t" + format + "\t" + sim + "\n")

