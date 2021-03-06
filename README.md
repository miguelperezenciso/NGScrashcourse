# NGScrashcourse
A few hour introductory course to NGS analysis

This is a few hour lecture that I teach at Universitat Autonoma of Barcelona [MSc in Bioinformatics](https://mscbioinformatics.uab.cat/). It contains first basic steps for someone who has not had exposure to NGS analyses, yet is familiar with linux commands. For simplicity, I base teaching in one of the smallest microorganism ([Mycoplasma genitalium](https://www.ncbi.nlm.nih.gov/genome/?term=Mycoplasma%20genitalium)) and a few simulated reads using [ART](https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm). 

### WARNINGGGGG
- **NGS data are massive yet noisy.** 
- **Caution and quality control are a must in every step, especially when analyzing several samples simultaneously.**

### Tools needed
Dozens of tools have been developed for NGS data, and each kind of data (eg, genomic, RNAseq, microbiome...) requires its own specific tools. It focuses on SNP calling. I may add new pipelines in the future for other type of analyses. 

Note: This is a samtools based pipeline, alternative pipelines are based on GATK. Hopefully, GATK pipelines will be added soon.

You need to install:
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/): This is a java tool that provides graphics of sequence quality. NGS consists of data in compressed [fastq](https://en.wikipedia.org/wiki/FASTQ_format) format.

- [BWA](https://sourceforge.net/projects/bio-bwa/files/). This is one of the most widely used aligner for genomic data.

	`git clone https://github.com/lh3/bwa.git`
	`cd bwa; make` 

- samtools, bcftools and hstslib: http://www.htslib.org/download/. A series of tools for manipulating bam files and variant calling. Once extracting, for each of them type 'make' and move executables to bin or to add to path. 
       
- [picard](http://broadinstitute.github.io/picard/). this is a series of tools to manipulate bam files and reads. Download and move **picard.jar** file to bin folder or add to path.

- [vcftools](https://sourceforge.net/projects/vcftools/). Filters and extract info from vcf files, format conversion with plink. Download, do 'make' and move executable to bin or to add to path.

- [igv](http://software.broadinstitute.org/software/igv/). this is a visualizer of bam files.

- [bedtools](https://bedtools.readthedocs.io/en/latest/). Check also bedtools, you will need it at some point during your career.

### Folder scheme
I suggest to have separate folders to organize the different analysis steps, but this is very personal. In the exercise, we have a folder **assembly** with the assembly and all required indices, a folder **reads** with all sequence data, bam files are stored in **bamfiles** directoy and vcf files in **varfiles**. Optionally, a **bin** folder contains executables.

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


### Definitions

	# define variables with program names or install in your path
	fastqc=$DIRBIN/FastQC/fastqc
	bwa=$DIRBIN/bwa
	samtools=$DIRBIN/samtools
	bcftools=$DIRBIN/bcftools
	picard=$DIRBIN/picard.jar
	vcftools=$DIRBIN/vcftools

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

### Step 0: Index reference genome

	cd $DIRASSEMBLY
	$bwa index -a bwtsw $DIRASSEMBLY/$ASSEMBLY
	
	# creating index with samtools
	$samtools faidx $ASSEMBLY

	# dictionary with picard tools' CreateSequenceDictionary (same name -> dict=reference)
	java -jar $picard CreateSequenceDictionary R=$DIRASSEMBLY/$ASSEMBLY O=$DIRASSEMBLY/$ASSEMBLY.dict
	# back
	cd $DWD

### Step 1: Check reads quality

	cd $DIRDATA
	zcat $DIRDATA/$READS_PE1 | $fastqc stdin
	zcat $DIRDATA/$READS_PE1 | $fastqc stdin
	# move back to wkng folder
	cd $DWD
	
### Step 2 BWA alignment and refinment
This is the most expensive part, both in CPU and memory usage

	# make a directory for that sample
	mkdir $DIRBAM/$OUT
	cd $DIRBAM/$OUT

	# This adds a tag to the alignment
	TAG="@RG\tID:$OUT\tSM:$OUT"
	time $bwa mem -t $NP -R $TAG $DIRASSEMBLY/$ASSEMBLY $DIRDATA/$READS_PE1 $DIRDATA/$READS_PE2 $samtools view -b - > $OUT.bam

	# sort
	time $samtools sort -O bam -T tmp $OUT.bam > $OUT.sort.bam
	time $samtools index $OUT.sort.bam

	# rm duplicates with picard, this removes potential PCR duplicates
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

### Step 3. samtools SNP CALLING with gVCF blocks 
This is the trickiest step, the one with most options and highly dependent on depth and sequence quality. Filtering SNPs is a must. Distrust indels more than SNPs. Beware of options with snp calling using several samples together. 

	cd $DIRVCF

	# check meaning of options typing "$bcftools"
	$bcftools mpileup -Ov -g $MINCOV -q $MAPQ -Q $BASEQ -d $MAXCOV -f $DIRASSEMBLY/$ASSEMBLY \
      		$DIRBAM/$OUT/$OUT.rmdup.bam | \
      		$bcftools call -g$MINCOV -mOz -o $OUT.gvcf.gz

	# filtering (see https://github.com/samtools/bcftools/wiki/HOWTOs#variant-filtering)
	$bcftools filter -O v -g3 -s LOWQUAL -e"%QUAL<$SNPQ || %MAX(INFO/DP)<$MINCOV || %MAX(INFO/DP)>$MAXCOV" \
      		$OUT.gvcf.gz | $bcftools view -f PASS -O z > $OUT.flt.gvcf.gz

	#--> EXERCISES
	#    Check the meaning of each option (eg, type $bcftools)
	#    Inspect the vcf file, what is each field?
	#    Do you find any indel?
	#    How many variants were filtered out (compare $OUT.flt.gvcf and $OUT.gvcf)? why?
	#    Count how many heterozygous snps, homozygous snps, why you find only one class?

	cd $DWD

### Step 4. Visualize some SNPs with IGV (http://software.broadinstitute.org/software/igv/) 

	# Start IGV (you need java8)
	# You need to load the Micoplasma genome in $DIRASSEMBLY
	java -Xmx1000m -jar $igv

### Bonus: Downloading sequences from SRA archive (https://www.ncbi.nlm.nih.gov/sra)
This is a most useful resource for NGS data. here I show how to automtaically download NGS reads from given samples. You need the SRR id.
**WARNING: this can take a lot of time and resources !!!!**

	# You need faster-qdump and aspera
	# To install aspera
	#    - https://www.ncbi.nlm.nih.gov/books/NBK242625/
	#    - http://downloads.asperasoft.com/connect2/

	# Aspera is a fast connection downloader    
	ASPERA=~/.aspera

	cd $DIRDATA
	# Choose a read to download, should start with SRR, SRX, ERR
	SRR=ERR4868557 # corresponds to a M genitalium sequenced with MiSeq
	# Exercise: inspect info about SRR6650027
	# this is the directory holding the compressed sequences
	DIRSRR=/sra/sra-instant/reads/ByRun/sra/${SRR:0:3}/${SRR:0:6}/$SRR
	# this downloads the sequences in $DIRDATA/SRR directory
	$ASPERA/connect/bin/ascp -i  $ASPERA/connect/etc/asperaweb_id_dsa.openssh \
                                 -k1 -Tr -l100m anonftp@ftp-private.ncbi.nlm.nih.gov:$DIRSRR $DIRDATA
	cd $SRR
	# uncompress into fastq, faster-qdump can be found in https://github.com/ncbi/sra-tools
	# this is actually the slowest step
	fasterq-dump -e $NP --split-files $SRR.sra -O $DIRDATA/$SRR
	rm $SRA.sra
	
	cd $DWD


### Exercises
 * Choose a different strain from experiment ERP004545. Note: you need to install aspera. Optionally, choose some of the already available five samples.
 * Perform alignment as described.
 * Perform multiple snp calling with about 5-10 samples together. 
 * Using vcftools, do a plot of allele frequency and missingness across all SNPs.
 * Identify common SNPs between samples.

### Exploring...
 - Get acquainted with major sequenicng technologies: Illumina (https://en.wikipedia.org/wiki/Illumina_dye_sequencing), Oxford Nanopore (https://en.wikipedia.org/wiki/Nanopore_sequencing), PacBio(https://en.wikipedia.org/wiki/Single-molecule_real-time_sequencing)...

 - Browse [GTAK](https://software.broadinstitute.org/gatk) best practices (https://software.broadinstitute.org/gatk/best-practices/):In particular, you can take a look at  https://software.broadinstitute.org/gatk/best-practices/workflow?id=11165; and https://software.broadinstitute.org/gatk/best-practices/workflow?id=11145
    
 - Browse the SRA archive (https://www.ncbi.nlm.nih.gov/sra/docs), which holds NGS data of numerous types and species. Project ERP004545 corresponds to M genitalium.
 
  - Nature Review Genetics series on [Applications of next-generation sequencing](https://www.nature.com/collections/jmgqdxpvsk)
  
  - ...
