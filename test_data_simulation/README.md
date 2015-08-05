Simulation Overview
===================

This file describes how to simulate an Alu polymorphism data set using the reference genome as input and writing bam files aligned with BWA mem as output.

You have two options:

1. Use the precomputed files `deletions.vcf` and `insertions.vcf` to re-generate the test data with 110 Alu polymorphisms on chr21 described in Qian et al. (2015). In this case you can skip step 1 below.
   Use hg19 in step 2 and the parameters as specified in the examples of step 3.
2. Generate your own set of insertions or deletions using the scripts `simulate_deletions.py` and `simulated_insertions.py`.
   Here, you need to provide a file with Alu positions in the reference genome.

Step 1 is only necessary if you choose option 1.
In both cases, you need to prepare the reference genome (step 2) and run the skript `simulate_bam.py` (step 3) for generating and mapping reads using MASON and BWA-mem.


Step 1 - Create a set of insertions/deletions
---------------------------------------------

The two python scripts `simulate_deletions.py` and `simulate_insertions.py` read the sequence of a reference chromosome from a fasta file as well as a file with Alu positions and create a VCF file and, in the case of insertions, a modified reference chromosome:
In addition, both scripts write a fasta file with the sequences of the selected Alu elements. These fasta files are not necessary for the following steps and you are free to ignore or delete them.

As input, the scripts expect only a single chromosome in the fasta file.
The file with Alu positions has to have six columns, but only columns 1-4 and 6 are being used:

    CHR  BEGIN  END  FAMILY  TYPE  ORIENTATION

For example:

    21  9986723   9987040   AluJb   SINE  +
    21  16253851  16254014  AluJo   SINE  +
    21  16347146  16347430  AluJb   SINE  +
    21  16735634  16735945  AluSx1  SINE  -
    21  17369195  17369499  AluY    SINE  +
    ...

If you run the scripts without parameters, they should print a short usage line:

    $ simulate_deletions.py
    Usage: ./simulate_deletions.py <ref.fa> <alu positions> [<num del> [<seed>]]
    
    $ simulate_insertions.py
    Usage: ./simulate_insertions.py <ref.fa> <alu positions> [<num ins> [<seed>]]

Replace `<ref.fa>` by your fasta file containing the chromosome and `<alu positions>` by the file listing the Alu positions.
Optionally, you can specify the number of insertions/deletions that you want to have in your data set, default is 100.
Further, you can specify a seed for the random number generator to be able to reproduce your data set.

Both python scripts output a VCF file with a set of Alu insertions or deletions.
In addition, `simulate_insertions.py` outputs a fasta file with a modified chromosome sequence.


Step 2 - Prepare the reference genome and BWA index
---------------------------------------------------

For creating a deletion data set, we only need a bwa index of the reference genome generated as usual:

    $ bwa index /path/to/genome.fa

For creating the insertion data set, you need to replace the modified chromosome in the reference genome before creating the bwa index.
You can, for example, use the following commands to do this:

* Determine the lines where the chromosome you want to replace starts and ends:


    $ grep -n ">" /path/to/genome.fa.fai
    ...
    39737155:>21
    40424726:>22
    ...

* Concatenate the reference genome with the head, cat and tail commands:


    $ head -n 39737154 /path/to/genome.fa > genome_ins.fa
    $ cat chr21.fa >> genome_ins.fa
    $ tail -n +40424726 /path/to/genome.fa >> genome_ins.fa

Make sure to substract 1 from the start line number in the head command.

Finally, create the bwa index for your modified reference genome:

    $ bwa index genome_ins.fa


Step 3 - Simulate and align the reads
-------------------------------------

The script `simulate_bam.sh` simulates the short read data for a single individual.
First, it creates two haplotypes with a subset of the Alu polymorphisms from the input file.
Then, it simulates the reads.
Finally, it alings the reads back to the reference genome using BWA mem.

You need to have installed BWA (e.g. version 0.7.10-r789), samtools (e.g. version 1.1), and Mason (e.g. version 2.0alpha1).
Before running the script, open the script in an editor and modify the `PATHS` section to match the locations of binaries on your system:

* `BWADIR` - Path to BWA mem, available at https://github.com/lh3/bwa
* `SAMDIR` - Path to samtools, available at https://github.com/samtools/samtools
* `SCRIPTDIR` - Path to add_variants.py, one of the PopAlu data simulation scripts
* `MASONDIR` - Path to the SeqAn tools mason_variator and mason_simulator available at http://www.seqan.de/projects/mason/

Before running the script, have a look at its five parameters described at the top of the script:

    WD=$1      # path to output directory for simulated individual (created if it does not exist)
    REFCHR=$2  # fasta file with the (modified) reference chromosome
    GENOME=$3  # fasta file with the (modified) whole reference genome (BWA index files have to be present in the same directory)
    INS=$4     # VCF file with insertions (deleted from reference) or deletions
    SEED=$5    # seed for simulation (use a different seed for each individual)

Run the script on the files created in the previous steps, e.g.

    $ simulate_bam.sh INS00 chr21_ins.fa genome_ins.fa insertions.vcf 0
    $ simulate_bam.sh INS01 chr21_ins.fa genome_ins.fa insertions.vcf 1
    ...
    $ simulate_bam.sh INS99 chr21_ins.fa genome_ins.fa insertions.vcf 99

for simulating an insertions data set or

    $ simulate_bam.sh DEL00 chr21.fa /path/to/genome.fa deletions.vcf 0
    $ simulate_bam.sh DEL01 chr21.fa /path/to/genome.fa deletions.vcf 1
    ...
    $ simulate_bam.sh DEL99 chr21.fa genome_ins.fa insertions.vcf 99

for simulating a deletions data set.

On success, the following files will be created in the specified output directories:

* `hap_0.vcf` and `hap_1.vcf`: VCF files with all variants inserted into the reference before simulating reads (including SNPs and indels in addition to a subset of the deletions/insertions from the input VCF file).
* `truth_0.bam` and `truth_1.bam`: The original locations (= correct alignments) of the reads to the reference.
* `mapped.bam`: An alignment of BWA of the simulated reads to the reference genome (= input to Alu polymorphism discovery programs).
