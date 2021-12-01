#!/bin/bash -x
#
# What is new: 2020-12-01 11:06:26 
#  - New GATK options for multisample calling
#  - GATK includes picard
#  - rm download from SRA, need fixing
#  - some notation changes

# What is new: 2018-11-30 11:01:46 
#  - option to download SRA sequences
#  - multisample variant calling, remove depth restriction
#  - include vcftools
#  - new exercises
# NGS bwa + samtools + GATK pipeline
# M Perez-Enciso
# miguel.perez@uab.es


# Exercise: browse the corresponding websites and get an idea of what they are for and general options
#     https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
#     http://bio-bwa.sourceforge.net/
#     https://software.broadinstitute.org/gatk/
#     http://software.broadinstitute.org/software/igv/
#     https://vcftools.github.io/index.html
#     https://bedtools.readthedocs.io/en/latest/


##############################################################################################
#                              EXECUTABLES                                                   #
#  download and install in bin directory
##############################################################################################
# FastQC download from https://www.bioinformatics.babraham.ac.uk/projects/download.html#fastqc

# bwa: https://sourceforge.net/projects/bio-bwa/files/
	git clone https://github.com/lh3/bwa.git
	cd bwa; make
# samtools, bcftools and hstslib: http://www.htslib.org/download/
# once extracting, for each of them 
        make
# and move executable to $DIRBIN
# GATK: download from https://github.com/broadinstitute/gatk/releases
# unzip and move to $DIRBIN
# vcftools: https://sourceforge.net/projects/vcftools/
        make
# and move executable to $DIRBIN
# igv: http://software.broadinstitute.org/software/igv/

Here it runs via web: https://igv.org/app/

##############################################################################################
#                              DIRECTORIES                                                   #
# this arrangement is simply a proposal, you can arrange as fits you best                    #
##############################################################################################
DWD=$(pwd)
# this contains the assembly
DIRASSEMBLY=$DWD/assembly
## Downloaded from https://www.ncbi.nlm.nih.gov/genome/474?genome_assembly_id=300158
ASSEMBLY=GCF_000027325.1_ASM2732v1_genomic
# this contains the reads
DIRDATA=$DWD/reads
# this contains the binaries
DIRBIN=$DWD/bin
# this will contain the alignment files
DIRBAM=$DWD/bamfiles
# this will contain the vcf files
DIRVCF=$DWD/varfiles


# define variables with program names or install in your path
fastqc=$DIRBIN/FastQC/fastqc
bwa=$DIRBIN/bwa
samtools=$DIRBIN/samtools
bcftools=$DIRBIN/bcftools
# adjust version as needed
gatk=$DIRBIN/gatk-4.2.3.0/gatk
vcftools=$DIRBIN/vcftools

# check what a variable contains
echo $vcftools

##############################################################################################
#                              SAMPLE NAME AND MAIN PARAMETERS                               #
##############################################################################################

################# global variables
MINCOV=5        # min coverage
MAXCOV=30       # maximum coverage computed later as 2xmean_depth+1
SNPQ=20         # min snp quality
MAPQ=20         # min map quality
BASEQ=20        # min base quality
NP=2            # no. of threads

# samples to be analyzed
sample1=Individual1
sample2=Individual2


################################################################################
# 0. INDEXING THE GENOME Build index sequence archive for reference sequence
################################################################################
   cd $DIRASSEMBLY
   $bwa index -a bwtsw $DIRASSEMBLY/$ASSEMBLY.fna

   # creating index with samtools
   $samtools faidx $ASSEMBLY.fna

   # dictionary with picard tools' CreateSequenceDictionary (same name -> dict=reference)
   $gatk CreateSequenceDictionary R=$DIRASSEMBLY/$ASSEMBLY.fna O=$DIRASSEMBLY/$ASSEMBLY.dict

   # back to working directory
   cd $DWD

################################################################################
# 1. CHECK QUALITY
################################################################################

  cd $DIRDATA
  # run interactively 
  $fastqc
  # or in batch , results should be in $DIRDATA folder
  zcat $DIRDATA/${sample1}.out.fas.1.fq.gz | $fastqc stdin
  zcat $DIRDATA/${sample1}.out.fas.2.fq.gz | $fastqc stdin
  cd $DWD

################################################################################
# 2. BWA ALIGNMENT AND GATK REALIGNMENT 
################################################################################
# reads are in files reads are in file ${sample}.out.fas.1.fq.gz and ${sample}.out.fas.2.fq.gz

   # make a directory for that sample
   mkdir $DIRBAM/$sample1
   cd $DIRBAM/$sample1

   # This adds a tag to the alignment
   TAG="@RG\tID:$sample1\tSM:$sample1"

   # Names of files containing paired end reads for individual $sample1
   READS_PE1=$DIRDATA/${sample1}.out.fas.1.fq.gz
   READS_PE2=$DIRDATA/${sample1}.out.fas.2.fq.gz

   # alignment, reads are in file ${sample1}.out.fas.1.fq.gz and ${sample1}.out.fas.2.fq.gz
   time $bwa mem -t $NP -R $TAG $DIRASSEMBLY/$ASSEMBLY.fna $READS_PE1 $READS_PE2 | \
        $samtools view -b - > $sample1.bam

   # sort
   time $samtools sort -O bam -T tmp $sample1.bam > $sample1.sort.bam
   time $samtools index $sample1.sort.bam

   # rm duplicates
   time $gatk MarkDuplicatesSpark \
                         -I $sample1.sort.bam \
                         -O $sample1.rmdup.bam \
                         --remove-sequencing-duplicates

   # recalibrate base quality with GATK if you have a list of known SNPs.
   # check https://software.broadinstitute.org/gatk/best-practices/workflow?id=11165

   # depth file
   time $samtools depth -q $BASEQ -Q $MAPQ $sample1.rmdup.bam | awk '{print $3}'  | \
                 sort | uniq -c | sort -n -k2 > $sample1.rmdup.depth
   # Exercise: find out what this instruction does
   # Exercise: what is the format of *depth file?

   # computes mean and max depth to be used, min depth is defined in $MINCOV
   Q=`awk '$2>0 {a+=$1*$2;n+=$1} END {print a/n}' "$DIRBAM/$sample1/$sample1.rmdup.depth" `
   # recommended maximum depth is twice average depth
   MAXCOV=`echo "tmp=2*($Q+0.5); tmp /= 1; tmp" | bc`
   echo 'sample meanDepth maxDepth ' $sample1 $Q $MAXCOV
   # Exercise: find out what instructions above do
   #           what is the bc shell instruction?

   #--> EXERCISE
   #    How many reads in each bam file:  samtools flagstat
   #    How many reads with map quality > 20:  samtools view -q 20 $OUT.bam | wc -l
   #    Perform alignment for sample Individual2 (sample2=Individual2)

   cd $DWD

################################################################################
# 3a. samtools SNP CALLING with gVCF blocks
################################################################################

   cd $DIRVCF

   # check meaning of options typing "$bcftools"
   $bcftools mpileup -Ov -g $MINCOV -q $MAPQ -Q $BASEQ -d $MAXCOV -f $DIRASSEMBLY/$ASSEMBLY.fna \
      $DIRBAM/$sample1/$sample1.rmdup.bam | \
      $bcftools call -g$MINCOV -mOz -o $sample1.gvcf.gz

   # filtering (see https://github.com/samtools/bcftools/wiki/HOWTOs#variant-filtering)
   $bcftools filter -O v -g3 -s LOWQUAL -e"%QUAL<$SNPQ || %MAX(INFO/DP)<$MINCOV || %MAX(INFO/DP)>$MAXCOV" \
      $sample1.gvcf.gz | $bcftools view -f PASS -O z > $sample1.flt.gvcf.gz

   # Multiple sample caling with samtools is performed adding multiple bam files in bcftools
   # either listing all bam files or specifying a bam list file with -b option
   # http://www.htslib.org/workflow/#mapping_to_variant
   # http://samtools.github.io/bcftools/bcftools.html#mpileup

   #--> EXERCISE
   #    Check the meaning of each option (eg, type $bcftools)
   #    Inspect the vcf file, what is each field?
   #    Do you find any indel?
   #    How many variants were filtered out (compare $sample1.flt.gvcf and $sample1.gvcf)? why?
   #    Count how many heterozygous snps, homozygous snps, why you find only one class?

   cd $DWD

################################################################################
# 3b. GATK Multiple sample SNP CALLING with gVCF blocks
################################################################################
# https://gatk.broadinstitute.org/hc/en-us/articles/360035890411?id=3893
# single sample calling

   cd $DIRVCF

   $gatk --java-options "-Xmx4g" HaplotypeCaller \
      -R $DIRASSEMBLY/$ASSEMBLY.fna \
      -I $DIRBAM/$sample1/$sample1.rmdup.bam \
      -O $sample1.g.vcf.gz \
      -ERC GVCF

   # You need to align reds for sample two
   $gatk --java-options "-Xmx4g" HaplotypeCaller \
      -R $DIRASSEMBLY/$ASSEMBLY.fna \
      -I $DIRBAM/$sample2/$sample2.rmdup.bam \
      -O $sample2.g.vcf.gz \
      -ERC GVCF

   # Combine gvcfs
   $gatk CombineGVCFs \
      -R $DIRASSEMBLY/$ASSEMBLY.fna \
      --variant $sample1.g.vcf.gz \
      --variant $sample2.g.vcf.gz \
      -O cohort.g.vcf.gz

   # Multisample calling
   $gatk --java-options "-Xmx4g" GenotypeGVCFs \
      -R $DIRASSEMBLY/$ASSEMBLY.fna \
      -V cohort.g.vcf.gz \
      -O cohort.vcf.gz

    cd $DWD

####################################################################################
# 4. Visualize some SNPs with IGV (http://software.broadinstitute.org/software/igv/) 
####################################################################################

# Start IGV (you need java8)
# You need to load the Micoplasma genome in $DIRASSEMBLY
java -Xmx1000m -jar $igv

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
#   4- Browse MiSeq, Oxford nanopore technologies, PacBio
####################################################################################





