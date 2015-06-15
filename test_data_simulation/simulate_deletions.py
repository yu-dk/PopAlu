#!/usr/bin/python

import random
from sys import argv

# -----------------------------------------------------------------------------
# initialize parameters
# -----------------------------------------------------------------------------

if len(argv) != 3 and len(argv) != 4 and len(argv) != 5:
   print "Usage:", argv[0], "<ref.fa> <alu positions> [<num del> [<seed>]]"
   exit()

# -----------------------------------------------------------------------------

chrfile = argv[1]
alufile = argv[2]

num_deletions = 100
if len(argv) > 3:
    num_deletions = int(argv[3])
    
seed = 0
if len(argv) > 4:
    seed = int(argv[4])

min_alu_len = 150
max_alu_len = 350



# -----------------------------------------------------------------------------
# functions
# -----------------------------------------------------------------------------

# write_seq()
#
# write the sequence <seq> to file <file> with line breaks every <line_width>
# characters

def write_seq(file, seq, line_width):
    for i in range(0, len(seq)/line_width + 1):
        pos = i*line_width
        file.write(seq[pos:pos+line_width] + "\n")

# -----------------------------------------------------------------------------



random.seed(seed)


# -----------------------------------------------------------------------------
# read the chromosome sequence
# -----------------------------------------------------------------------------

chr_id = ""
chr_seq = ""
for line in open(chrfile):
    if line[0] == ">" :
        if chr_id != "":
            print "ERROR: Expecting only one fasta record in input (" + chrfile + ")."
            exit()
        chr_id = line[1:-1]
        continue
    chr_seq += line[:-1]

print "Length of original reference: " + str(len(chr_seq))
chr_seq = chr_seq.upper()


# -----------------------------------------------------------------------------
# read the alu positions and pick num_deletions of them
# -----------------------------------------------------------------------------

alu_pos = []
for line in open(alufile):
    if line[0:3] != "chr":
        print "ERROR_ Expecting lines in alufile to start with 'chr'."
        exit()

    # read the begin and end position from line
    line = line.split()
    begin = int(line[1])
    end = int(line[2])
    if end-begin >= min_alu_len and end-begin <= max_alu_len:
        alu_pos.append([begin, end, line[3], line[5]])

num_alu = len(alu_pos)
positions = []

for i in range(0, num_deletions):

    # pick random position on chr
    pos = random.randint(0, num_alu)

    positions.append(alu_pos[pos])

# sort the position in ascending order
positions.sort()

# -----------------------------------------------------------------------------
# write segments as deletions
# -----------------------------------------------------------------------------

vcffile = open("deletions.vcf", "w")
fafile = open("deletions.fa", "w")

# write header
vcffile.write("##fileformat=VCFv4.1\n")
vcffile.write("##source=" + argv[0] + "\n")
vcffile.write("##reference=" + chrfile + "\n")
vcffile.write("##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">\n")
vcffile.write("##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">\n")
vcffile.write("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n")
vcffile.write("##contig=<ID=" + chr_id + ",length=" + str(len(chr_seq)) + ">\n")
vcffile.write("#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  simulated)\n")

# write records
for i in range(0, num_deletions):
    begin = positions[i][0]
    end = positions[i][1]
    type  = positions[i][2]
    orientation  = positions[i][3]
    
    # pick random frequency
    freq = random.random()

    # substring from begin to end position
    ins_seq = chr_seq[begin:end]
    
    vcffile.write(chr_id + "\t")
    vcffile.write(str(begin) + "\t")
    vcffile.write(type + "_chr21_" + str(begin) + "\t")
    vcffile.write(chr_seq[begin-1] + "\t")
    vcffile.write(chr_seq[begin-1] + ins_seq + "\t")
    vcffile.write("." + "\t")
    vcffile.write("PASS" + "\t")
    vcffile.write("SVTYPE=DEL;SVLEN=-" + str(len(ins_seq)) + ";")
    vcffile.write("AF=" + "{0:.4f}".format(freq) + "\t")
    vcffile.write("." + "\t")
    vcffile.write("1" + "\n")
    
    fafile.write(">chr:" + chr_id + "|")
    fafile.write("pos:" + str(begin) + " ")
    fafile.write("len:" + str(len(ins_seq)) + " ")
    fafile.write("ori:" + str(orientation) + " ")
    fafile.write("type:" + str(type) + " ")
    fafile.write("freq:" + "{0:.4f}".format(freq) + "\n")
    write_seq(fafile, ins_seq, 70)

vcffile.close()
fafile.close()

