#!/bin/bash

###############################################################################
# INPUT - PROGRAM PARAMETERS

WD=$1      # path to directory for simulated individual (will be created)
REFCHR=$2  # fasta file with the (modified) reference chromosome
GENOME=$3  # fasta file with the (modified) whole reference genome (BWA index
           #    files have to be present in the same directory)
INS=$4     # VCF file with insertions (deleted from reference) or deletions
SEED=$5    # seed for simulation (use a different seed for each individual)

# IMPORTANT:
#  - Use the modified reference chromosome created by simulate_inserts.py for
#    simulating data with insertions.
#  - Use the original reference chromosome and genome for simulating data with
#    deletions.

COV=25    # the read coverage

###############################################################################
# OUTPUT

# If this script finishes successfully, WD will contain
#   - hap_0.fa, hap_1.fa
#   - truth_0.bam, truth_1.bam
#   - mapped.bam, mapped.bam.bai


###############################################################################
# PATHS

BWADIR=~/bin                                        # Path to bwa
SAMDIR=~/bin                                        # Path to samtools
SCRIPTDIR=~/tools/PopAlu/test_data_simulation       # Path to add_variants.py
MASONDIR=~/bin                                      # Path to mason


###############################################################################
# trap error handler: print location of last error and exit

function my_trap_handler()
{
        MYSELF="$0"              # equals to my script name
        LASTLINE="$1"            # argument 1: last line of error occurence
        LASTERR="$2"             # argument 2: error code of last command
        echo "${MYSELF}: line ${LASTLINE}: exit status of last command: ${LASTERR}"
        exit
}

# trap commands with non-zero exit code
trap 'my_trap_handler ${LINENO} $?' ERR


###############################################################################
# CREATE WORKING DIRECTORY

echo $(date) "Working directory is ${WD}"

if [ ! -d ${WD} ]; then
    mkdir -p ${WD}
fi


###############################################################################
# CREATE TWO HAPLOTYPES

echo $(date) "Creating haplotypes"

# different seeds for the two haplotypes
SEED1=$(( ${SEED} * 2 ))
SEED2=$(( ${SEED1} + 1 ))

# simulate SNPs and small indels
${MASONDIR}/mason_variator --seed ${SEED1} --in-reference ${REFCHR} --out-vcf ${WD}/hap_0.vcf --sv-indel-rate 0 --sv-inversion-rate 0 --sv-translocation-rate 0 --sv-duplication-rate 0
${MASONDIR}/mason_variator --seed ${SEED2} --in-reference ${REFCHR} --out-vcf ${WD}/hap_1.vcf --sv-indel-rate 0 --sv-inversion-rate 0 --sv-translocation-rate 0 --sv-duplication-rate 0

# add insertions to vcf file of SNPs and indels
${SCRIPTDIR}/add_variants.py ${INS} ${WD}/hap_0.vcf ${SEED1}
${SCRIPTDIR}/add_variants.py ${INS} ${WD}/hap_1.vcf ${SEED2}

# sort the vcf files by position (without affecting the header)
(grep "^#" ${WD}/hap_0.vcf && grep -v "^#" ${WD}/hap_0.vcf | sort -k2,2n) > ${WD}/hap_0_sorted.vcf
(grep "^#" ${WD}/hap_1.vcf && grep -v "^#" ${WD}/hap_1.vcf | sort -k2,2n) > ${WD}/hap_1_sorted.vcf
mv ${WD}/hap_0_sorted.vcf ${WD}/hap_0.vcf
mv ${WD}/hap_1_sorted.vcf ${WD}/hap_1.vcf


###############################################################################
# SIMULATE READS

echo $(date) "Simulating reads"

CHARS=`grep -v ">" ${REFCHR} | wc -c`
LINES=`grep -v ">" ${REFCHR} | wc -l`
LEN=$(( ${CHARS} - ${LINES} ))
NUM_FRAGS=$(( ${LEN}/101*${COV}/4 )) # chromosome length / read length * coverage / 4
                                     # divide by 4 since we have two haplotypes and two reads per fragment

# reads for first haplotype
${MASONDIR}/mason_simulator --seed ${SEED1} --num-fragments ${NUM_FRAGS} --out ${WD}/reads.0.fastq --out-right ${WD}/reads.1.fastq --out-alignment ${WD}/truth_0.bam --input-reference ${REFCHR} --input-vcf ${WD}/hap_0.vcf --fragment-mean-size 500 --read-name-prefix "hap:${SEED1}|sim:" --illumina-read-length 101 --illumina-prob-insert 0.0001 --illumina-prob-deletion 0.0001

# reads for second haplotype in 10 batches
${MASONDIR}/mason_simulator --seed ${SEED2} --num-fragments ${NUM_FRAGS} --out ${WD}/reads_1.0.fastq --out-right ${WD}/reads_1.1.fastq --out-alignment ${WD}/truth_1.bam --input-reference ${REFCHR} --input-vcf ${WD}/hap_1.vcf --fragment-mean-size 500 --read-name-prefix "hap:${SEED2}|sim:" --illumina-read-length 101 --illumina-prob-insert 0.0001 --illumina-prob-deletion 0.0001

# append reads of second haplotype to reads of first haplotype
cat ${WD}/reads_1.0.fastq >> ${WD}/reads.0.fastq
rm ${WD}/reads_1.0.fastq
cat ${WD}/reads_1.1.fastq >> ${WD}/reads.1.fastq
rm ${WD}/reads_1.1.fastq


###############################################################################
# MAP SIMULATED READS BACK TO WHOLE GENOME

echo $(date) "Mapping reads"

${BWADIR}/bwa mem ${GENOME} ${WD}/reads.0.fastq ${WD}/reads.1.fastq > ${WD}/mapped.sam
rm ${WD}/reads.0.fastq ${WD}/reads.1.fastq

# convert to sam, sort and index
${SAMDIR}/samtools view ${WD}/mapped.sam -S -h -b > ${WD}/mapped_unsorted.bam
rm ${WD}/mapped.sam
${SAMDIR}/samtools sort ${WD}/mapped_unsorted.bam ${WD}/mapped
rm ${WD}/mapped_unsorted.bam
${SAMDIR}/samtools index ${WD}/mapped.bam


###############################################################################

echo $(date) "Finished ${WD} successfully."

###############################################################################

