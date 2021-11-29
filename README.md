# NGS Crash Course
A few hour introductory course to NGS analysis, mainly for variant calling.

This is a few hour lecture that I usually teach at [Universitat Autonoma of Barcelona](www.uab.es) [MSc in Bioinformatics](https://mscbioinformatics.uab.cat/). It contains first basic steps for someone who has not had exposure to NGS analyses, yet is familiar with linux commands. For simplicity, I base teaching in one of the smallest microorganism ([Mycoplasma genitalium](https://www.ncbi.nlm.nih.gov/genome/?term=Mycoplasma%20genitalium)) and a few simulated reads using [ART](https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm). 

There exists dozens of procedures specific to many of the applications of NGS data.

**These notes are intended for variant calling with genome data.**

### WARNINGGGGG 1

This is a rapidly evolving field. 

### Why NGS? 

* To be honest, often because it is *very cheap*.
* Because one technology fits all:
	- De novo sequencing: assemblies of new species.
	- Re sequencing: variant discovery, both SNPs and structural variants
	- mRNA resequencing: discovery of new splicing events, quantification of transcription activity
	- microRNA characterization
	- Bisulfide sequencing: Epigenomics 
	- CHip technology: DNA motifs
	- Metagenomics
	- Medical applications: cancer genomes …
	- ...

Perhaps you can check the [encode project](https://www.encodeproject.org/) to find about the numerous applications of NGS.

### Do we really need NGS?
Perhaps we don't but is sooooo cheap that sequencing whatever has become irresistible!

Note sequence data is highly redundant and that you will never get many samples with complete sequence. NGS data are full of 'holes' (missing posiitons) even at very high coverage.

For many applications, high density genotyping is may be enough.

### WARNINGGGGG 2
- **NGS data are massive yet VERY NOISY.** 
- **You will be puzzled about how results can change depending on bioinformatics options utilized.**
- **Do NOT outsource this step. DO IT YOURSELF or someone you know well.**
- **Caution and quality control are a must in every step, especially when analyzing several samples simultaneously.**
- **Missing data are unavoidable.**

## QUALITY CONTROL IS A MUST IN EVERY NGS ANALYSIS STEP 
### VISUALIZATION TOOLS CAN REALLY HELP.

### What techniques are there?
Classical Sanger sequencing was a 'First' generation technology. Second or next generation sequencing was initially [454 technique](https://en.wikipedia.org/wiki/454_Life_Sciences) quickly supersed by [Illumina's technologies](https://en.wikipedia.org/wiki/Illumina_dye_sequencing). More recently, the third generation sequencing technologies deliver much longer reads in lower quantity and of lower quality. Representatives are [PacBio](https://www.sciencedirect.com/science/article/pii/S1672022915001345) and [Oxford nanopore](https://en.wikipedia.org/wiki/Nanopore_sequencing).

Each of these technologies may require slightly different software / options but bioinformatics pipelines are essentially the same.

### Which are the main steps?
Sequence data need to be aligned against a reference genome. Sequence quality check is warranted. Next, the alignment file needs to be sorted and can be polished by removing PCR duplicates and quality can be recalibrated with a known set of variants. 

There exist several options for variant calling, main ones are [samtools](http://www.htslib.org/download/) and [GATK](https://gatk.broadinstitute.org/hc/en-us) based.

![Figure 1](https://github.com/miguelperezenciso/NGScrashcourse/blob/master/steps.png)


### Main data formats
Fortunately, a set of standardized formats allow communication between softwares and pipelines. Main ones are **fastq** for sequence data, **BAM / CRAM** for aligned sequences and **vcf / gvcf** for polymorphisms. Many of these formats and softwares were developed for the [1000 human genome project](https://www.internationalgenome.org/). 

![](https://github.com/miguelperezenciso/NGScrashcourse/blob/master/formats.png)


**[fastq](https://en.wikipedia.org/wiki/FASTQ_format)**: It is a modified fasta with quality information, and it is the format that sequencing technologies. deliver their data.

	@SEQ_ID
	GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT	
	+
	!''*((((***+))%%%++)(%%%%).1***-+*''))**55CCF>>>>>>CCCCCCC65
	
- Line 1 begins with a '@' character and is followed by a sequence identifier and an optional description (like a FASTA title line).
- Line 2 is the raw sequence letters.
- Line 3 begins with a '+' character and is optionally followed by the same sequence identifier (and any description) again.
- Line 4 encodes the quality values for the sequence in Line 2, and must contain the same number of symbols as letters in the sequence.

Quality values (Q) are encoded using the so called 'Phred' scores (P) using [ASCII](https://es.wikipedia.org/wiki/ASCII) codes. The beauty of this is that numbers with several digits can be represented with a single letter. For instance, symbol '!' has ASCII code 33. Illumina 1.8+ uses Phred+33, raw reads typically (0, 41). The relation between Q and P 

	Q = -10 log10(P)
	
	P = 10^(-Q/10)

In python, function **ord** prints ASCII code, so P of error coded with symbol 'X' would be

	P = 10**(-(ord("X")-33)/10))
	
### EXERCISE 1:
	Write a python script that returns average P from a fastq file.
	
**BAM** format is the compressed version of [SAM](https://en.wikipedia.org/wiki/SAM_(file_format)) or Sequence Alignment Map. **SAM** format contains information on quality and position alignment for every read. You can visualize a BAM file with
	
	# visualize first lines of a BAM file (you do not normally need to do this)
	samtools view file.bam | head
	
**[vcf](https://samtools.github.io/hts-specs/VCFv4.2.pdf)**: format contains variants found for one or more samples. It contains differences between the sample and the reference genome. Both SNPs and indels can be represented. Example:

	##fileformat=VCFv4.2
	##fileDate=20090805
	##source=myImputationProgramV3.1
	##reference=file:///seq/references/1000GenomesPilot-NCBI36.fasta
	##contig=<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species="Homo sapiens",taxonomy=x>
	##phasing=partial
	##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
	##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
	##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
	##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">
	##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership, build 129">
	##INFO=<ID=H2,Number=0,Type=Flag,Description="HapMap2 membership">
	##FILTER=<ID=q10,Description="Quality below 10">
	##FILTER=<ID=s50,Description="Less than 50% of samples have data">
	##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
	##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
	##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
	##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">
	#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT NA00001 NA00002 NA00003
	20 14370 rs6054257 G A 29 PASS NS=3;DP=14;AF=0.5;DB;H2 GT:GQ:DP:HQ 0|0:48:1:51,51 1|0:48:8:51,51 1/1:43:5:.,.
	20 17330 . T A 3 q10 NS=3;DP=11;AF=0.017 GT:GQ:DP:HQ 0|0:49:3:58,50 0|1:3:5:65,3 0/0:41:3
	20 1110696 rs6040355 A G,T 67 PASS NS=2;DP=10;AF=0.333,0.667;AA=T;DB GT:GQ:DP:HQ 1|2:21:6:23,27 2|1:2:0:18,2 2/2:35:4
	20 1230237 . T . 47 PASS NS=3;DP=13;AA=T GT:GQ:DP:HQ 0|0:54:7:56,60 0|0:48:4:51,51 0/0:61:2
	20 1234567 microsat1 GTC G,GTCT 50 PASS NS=3;DP=9;AA=G GT:GQ:DP 0/1:35:4 0/2:17:2 1/1:40:3

Explanation:
* First lines start with **##** and are comments
* Line starting with **#CHROM** contains field names: chromosome, position, variant id (not usually available), base in reference, alternative base in sample(s), variant quality phred, whether variant reliable or not, information (see ## fields), format (see ## fields), sample name(s).
* Genotype is specified in columns 10 and successive, one for each sample; `0/0` means homozygous genotype as in reference, `0/1` is heterozygous, `1/1` is homozygous for alternative allele; `0|1` means phased genotype, as opposed to default `0/1`.
* Numbers after genotype code mean probability of each genotype

**[gvcf](https://github.com/broadinstitute/gatk-docs/blob/master/gatk3-faqs/What_is_a_GVCF_and_how_is_it_different_from_a_'regular'_VCF%3F.md)** format is defined only for single sample genotypes and contains information on which genome tracts are identical between sample and reference genome. useful for multisample variant calling.

**[gff3](http://www.ensembl.org/info/website/upload/gff3.html)** generic format to host annotations.

**WARNING 3**: Separator is tab (\t) by default in all these formats.

### Main tools needed
Dozens of tools have been developed for NGS data, and each kind of data (eg, genomic, RNAseq, microbiome...) requires its own specific tools. Here I focus on SNP calling. I may add new pipelines in the future for other type of analyses. 

Note: This is a samtools and GATK based pipeline. GATK is better documented (perhaps too much!), but is slower and not necessarily resulting in lower false discovery rate and higher power. 

![](https://github.com/miguelperezenciso/NGScrashcourse/blob/master/software.png)

You need to install:
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/): This is a java tool that provides graphics of sequence quality. NGS consists of data in compressed [fastq](https://en.wikipedia.org/wiki/FASTQ_format) format.

- [BWA](https://sourceforge.net/projects/bio-bwa/files/). This is one of the most widely used aligner for genomic data.

	`git clone https://github.com/lh3/bwa.git`
	`cd bwa; make` 

- samtools, bcftools and hstslib: http://www.htslib.org/download/. A series of tools for manipulating bam files and variant calling. Once extracting, for each of them type 'make' and move executables to bin or to add to path. 

- [GATK](https://gatk.broadinstitute.org/hc/en-us).

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
