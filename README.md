# NGScrashcourse
A few hour introductory course to NGS analysis

This is a few hour lecture that I teach at Universitat Autonoma of Barcelona [MSc in Bioinformatics](https://mscbioinformatics.uab.cat/). It contains first basic steps for someone who has not had exposure to NGS analyses, yet is familiar with linux commands. For simplicity, I base teaching in one of the smallest microorganism ([Mycoplasma genitalium](https://www.ncbi.nlm.nih.gov/genome/?term=Mycoplasma%20genitalium)) and a few simulated reads using [ART](https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm). 

### WARNINGGGGG
- **NGS data are massive yet noisy.** 
- **Caution and quality control are a must in every step, especially when analyzing several samples simultaneously.**

### Tools needed
Hundreds of tools have been developed for NGS data, and each kind of data (eg, genomic, RNAseq, microbiome...) requires its own specific tools. It focuses on SNP calling. I may add new pipelines in the future for other type of analyses. 

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

### Definitions

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


### Step 0: Check reads quality 

### Exploring...
 - Get acquainted with major sequenicng technologies: Illumina (https://en.wikipedia.org/wiki/Illumina_dye_sequencing), Oxford Nanopore (https://en.wikipedia.org/wiki/Nanopore_sequencing), PacBio(https://en.wikipedia.org/wiki/Single-molecule_real-time_sequencing)...

 - Browse [GTAK](https://software.broadinstitute.org/gatk) best practices (https://software.broadinstitute.org/gatk/best-practices/):In particular, you can take a look at  https://software.broadinstitute.org/gatk/best-practices/workflow?id=11165; and https://software.broadinstitute.org/gatk/best-practices/workflow?id=11145
    
 - Browse the SRA archive (https://www.ncbi.nlm.nih.gov/sra/docs), which holds NGS data of numerous types and species. Project ERP004545 corresponds to M genitalium.

## Exercises
 * Choose a different strain from experiment ERP004545, you need to install aspera

 * Run fastqc, Perform alignment as described

 * Perform multiple snp calling with about 5-10 samples together

4 Using vcftools, do a plot of allele frequency and missingness across all SNPs
