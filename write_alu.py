"""
input is rmsk data, eg, http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz
output is alu_chr* files
"""
import sys, gzip

f_in = sys.argv[1]  ###'eg: rmsk.txt.gz'

chrns = ['chr%d'%i for i in range(1, 23)] + ['chrX', 'chrY']
chrns_info = dict((i, {}) for i in chrns)

fin = gzip.open(f_in, 'rb')
for line in fin:
    tmp = line.split()
    chrn, pos, pos2, alu, strand = tmp[5], tmp[6], tmp[7], tmp[12], tmp[9]
    if 'Alu' not in alu:
        continue
    if not chrn in chrns_info:
        continue
    chrns_info[chrn][int(pos)] = '\t'.join([chrn, pos, pos2, alu, 'na', strand])
#    print chrn, pos
fin.close()

path1 = '/'.join(f_in.split('/')[:-1])
for chrn in chrns:
    with open(path1  + '/alu_' + chrn, 'w') as fout:
        for k in sorted(chrns_info[chrn].keys()):
            print >>fout, chrns_info[chrn][k]

print 'input', f_in
print 'output to %s/alu_chr*'%path1

