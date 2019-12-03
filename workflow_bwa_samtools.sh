#!/bin/bash -x
# NGS bwa + samtools + picard pipeline
# M Perez-Enciso
# mperezenciso@gmail.com

##############################################################################################
#                              DIRECTORIES                                                   #
# this arrangement is simply a proposal, you can arrange as fits you best                    #
##############################################################################################
# current folder
DWD=$(pwd)
# this should contain the assembly
DIRASSEMBLY=$DWD/assembly
## Downloaded from https://www.ncbi.nlm.nih.gov/genome/474?genome_assembly_id=300158
ASSEMBLY=GCF_000027325.1_ASM2732v1_genomic.fna
# this should contain the reads
DIRDATA=$DWD/reads
# this contains the binaries, alternatively, they can be accessed via default path
DIRBIN=$DWD/bin
# this will contain the alignment files
DIRBAM=$DWD/bamfiles
mkdir $DIRBAM
# this will contain the vcf files
DIRVCF=$DWD/varfiles
mkdir $DIRVCF

# Exercise: browse the corresponding websites and get an idea of what they are for and general options
#     https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
#     http://bio-bwa.sourceforge.net/
#     http://www.htslib.org/
#     https://software.broadinstitute.org/gatk/
#     http://broadinstitute.github.io/picard/
#     http://software.broadinstitute.org/software/igv/
#     https://vcftools.github.io/index.html
# Check also bedtools, you will need it at some point during your career
#     https://bedtools.readthedocs.io/en/latest/


##############################################################################################
#                              EXECUTABLES                                                   #
#                   download and install in bin directory                                    #
##############################################################################################
# FastQC download from https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
# bwa: https://sourceforge.net/projects/bio-bwa/files/
	git clone https://github.com/lh3/bwa.git
	cd bwa; make
# samtools, bcftools and hstslib: http://www.htslib.org/download/
# once extracting, for each of them 
        make
# and move executable to $DIRBIN
# picard: http://broadinstitute.github.io/picard/
# and move picard.jar to $DIRBIN

# vcftools: https://sourceforge.net/projects/vcftools/
        make
# and move executable to $DIRBIN

# igv: http://software.broadinstitute.org/software/igv/

# define variables with program names or install in your path
fastqc=$DIRBIN/FastQC/fastqc
bwa=$DIRBIN/bwa
samtools=$DIRBIN/samtools
bcftools=$DIRBIN/bcftools
picard=$DIRBIN/picard.jar
vcftools=$DIRBIN/vcftools


##############################################################################################
#                              SAMPLE NAME AND MAIN PARAMETERS                               #
##############################################################################################

# sample to be analyzed
OUT=Individual1

# reads from sample, they are in $DIRDATA
READS_PE1=${OUT}.out.fas.1.fq.gz
READS_PE2=${OUT}.out.fas.2.fq.gz

# global variables
MINCOV=5        # min coverage
MAXCOV=30       # maximum coverage computed later as 2xmean_depth+1
SNPQ=10         # min snp quality
MAPQ=20         # min map quality
BASEQ=20        # min base quality
NP=2            # no. of threads


################################################################################
# 0. INDEXING THE GENOME Build index sequence archive for reference sequence
################################################################################
   cd $DIRASSEMBLY
   $bwa index -a bwtsw $DIRASSEMBLY/$ASSEMBLY

   # creating index with samtools
   samtools faidx $ASSEMBLY

   # dictionary with picard tools' CreateSequenceDictionary (same name -> dict=reference)
   java -jar $picard CreateSequenceDictionary R=$DIRASSEMBLY/$ASSEMBLY O=$DIRASSEMBLY/$ASSEMBLY.dict
   cd $DWD

################################################################################
# 1. CHECK QUALITY
################################################################################

  cd $DIRDATA
  # run interactively 
  $fastqc
  # or in batch 
  zcat $DIRDATA/$READS_PE1 | $fastqc stdin
  zcat $DIRDATA/$READS_PE1 | $fastqc stdin
  cd $DWD

################################################################################
# 2. BWA ALIGNMENT AND GATK REALIGNMENT (~5-12h)
################################################################################

   # make a directory for that sample
   mkdir $DIRBAM/$OUT
   cd $DIRBAM/$OUT

   # This adds a tag to the alignment
   TAG="@RG\tID:$OUT\tSM:$OUT"
   time $bwa mem -t $NP -R $TAG $DIRASSEMBLY/$ASSEMBLY $DIRDATA/$READS_PE1 $DIRDATA/$READS_PE2 | \
        $samtools view -b - > $OUT.bam

   # sort
   time $samtools sort -O bam -T tmp $OUT.bam > $OUT.sort.bam
   time $samtools index $OUT.sort.bam

   # rm duplicates with picard
   time java -jar $picard MarkDuplicates \
                         REMOVE_DUPLICATES=true \
                         INPUT=$OUT.sort.bam \
                         OUTPUT=$OUT.rmdup.bam \
                         METRICS_FILE=metrics.out

   # recalibrate base quality with GATK if you have a list of known SNPs.
   # check https://software.broadinstitute.org/gatk/best-practices/workflow?id=11165

   # depth file
   time $samtools depth -q $BASEQ -Q $MAPQ $OUT.rmdup.bam | awk '{print $3}'  | \
                 sort | uniq -c | sort -n -k2 > $OUT.rmdup.depth
   # Exercise: find out what this instruction does

   # new indexed file
   $samtools index $OUT.rmdup.bam

   # computes mean and max depth to be used, min depth is defined in $MINCOV
   Q=`awk '$2>0 {a+=$1*$2;n+=$1} END {print a/n}' "$DIRBAM/$OUT/$OUT.rmdup.depth" `
   # recommended maximum depth is twice average depth
   MAXCOV=`echo "tmp=2*($Q+0.5); tmp /= 1; tmp" | bc`
   echo 'sample meanDepth maxDepth ' $OUT $Q $MAXCOV

   #--> EXERCISE
   #    How many reads in each bam file:  samtools flagstat
   #    How many reads with map quality > 20:  samtools view -q 20 $OUT.bam | wc -l

   cd $DWD

################################################################################
# 3. $samtools SNP CALLING with gVCF blocks 
################################################################################

   cd $DIRVCF

   # check meaning of options typing "$bcftools"
   $bcftools mpileup -Ov -g $MINCOV -q $MAPQ -Q $BASEQ -d $MAXCOV -f $DIRASSEMBLY/$ASSEMBLY \
      $DIRBAM/$OUT/$OUT.rmdup.bam | \
      $bcftools call -g$MINCOV -mOz -o $OUT.gvcf.gz

   # filtering (see https://github.com/samtools/bcftools/wiki/HOWTOs#variant-filtering)
   $bcftools filter -O v -g3 -s LOWQUAL -e"%QUAL<$SNPQ || %MAX(INFO/DP)<$MINCOV || %MAX(INFO/DP)>$MAXCOV" \
      $OUT.gvcf.gz | $bcftools view -f PASS -O z > $OUT.flt.gvcf.gz

   #--> EXERCISE
   #    Check the meaning of each option (eg, type $bcftools)
   #    Inspect the vcf file, what is each field?
   #    Do you find any indel?
   #    How many variants were filtered out (compare $OUT.flt.gvcf and $OUT.gvcf)? why?
   #    Count how many heterozygous snps, homozygous snps, why you find only one class?

   cd $DWD

####################################################################################
# 4. Visualize some SNPs with IGV (http://software.broadinstitute.org/software/igv/) 
####################################################################################

# Start IGV (you need java8)
# You need to load the Micoplasma genome in $DIRASSEMBLY
java -Xmx1000m -jar $igv

########################################################################################
# 5. OPTIONAL: Downloading sequences from SRA archive (https://www.ncbi.nlm.nih.gov/sra) 
########################################################################################
    # WARNING: this can take a lot of time and resources !!!!
    # You need faster-qdump and aspera
    # To install aspera
    #    - https://www.ncbi.nlm.nih.gov/books/NBK242625/
    #    - http://downloads.asperasoft.com/connect2/

    
    ASPERA=~/.aspera

    cd $DIRDATA
    # Choose a read to download, should start with SRR, SRX, ERRR
    SRR=ERR4868557 # corresponds to a M genitalium sequenced with MiSeq
    # inspect info about SRR6650027
    # this is the directory holding the compressed sequences
    DIRSRR=/sra/sra-instant/reads/ByRun/sra/${SRR:0:3}/${SRR:0:6}/$SRR
    # this downloads the sequences in $DIRDATA/SRR directory
    $ASPERA/connect/bin/ascp -i  $ASPERA/connect/etc/asperaweb_id_dsa.openssh \
                                 -k1 -Tr -l100m anonftp@ftp-private.ncbi.nlm.nih.gov:$DIRSRR $DIRDATA
    cd $SRR
    # uncompress into fastq, faster-qdump can be found in https://github.com/ncbi/sra-tools
    fasterq-dump -e $NP --split-files $SRR.sra -O $DIRDATA/$SRR
    rm $SRA.sra

exit 0

####################################################################################
# Exploring
####################################################################################
#   1- Browse GTAK best practices (https://software.broadinstitute.org/gatk/best-practices/)
#      - https://software.broadinstitute.org/gatk/best-practices/workflow?id=11165
#      - https://software.broadinstitute.org/gatk/best-practices/workflow?id=11145
#   2- Browse SRA archive
#      - https://www.ncbi.nlm.nih.gov/sra/docs
#      - SRA ERP004545 corresponds to M genitalium
#   3- Browse galaxy platform
#      - https://usegalaxy.org/
#   4- Browse MiSeq, Oxford nanopore technologies
####################################################################################

####################################################################################
# Exercises
####################################################################################
#   1- Each of you choose a different strain from experiment ERP004545, you need to install aspera
#   2- Run fastqc, Perform alignment as described, each of you individually
#   3- get together in groups and perform multiple snp calling with about 5-10 samples together
#      Note: remove the maximum depth option     
#   4- Using vcftools, do a plot of allele frequency and missingness across all SNPs
####################################################################################




