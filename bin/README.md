# INSTALLING TOOLS

Dozens of tools have been developed for NGS data, and each kind of data (eg, genomic, RNAseq, microbiome...) requires its own specific tools. Here I focus on SNP calling. 

![](https://github.com/miguelperezenciso/NGScrashcourse/blob/master/software.png)

First, move to **$DIRBIN** folder where I install all tools for convenience.

You need to install:

[**FastQC**](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/): This is a java tool that provides graphics of sequence quality. 
Documentation can be found [here](https://raw.githubusercontent.com/s-andrews/FastQC/master/INSTALL.txt).
You can run FastQC in one of two modes, either as an interactive graphical application.

	# running interactively in linux
	./fastqc
	# and navigate to load a compressed fastq file

	# run as shell
	zcat *fastq.gz | fastqc stdin

[**BWA**](https://sourceforge.net/projects/bio-bwa/files/). This is one of the most widely used aligner for genomic data.

	git clone https://github.com/lh3/bwa.git
	cd bwa
	make 

[**samtools, bcftools and hstslib**](http://www.htslib.org/download/). A series of tools for manipulating bam files and variant calling. Once extracted, for each of them type 'make' and move executables to bin or add to path.

	# download from http://www.htslib.org/download
	cd samtools-1.x    # and similarly for bcftools and htslib
	./configure --prefix=$DIRBIN
	make
	make install

[**GATK**](https://gatk.broadinstitute.org/hc/en-us).

	# download from https://github.com/broadinstitute/gatk/releases
	# unzip and move to $DIRBIN

[**vcftools**](https://sourceforge.net/projects/vcftools/). Filters and extract info from vcf files, format conversion with plink. 

	# download from https://sourceforge.net/projects/vcftools
	make 
	# move executable to bin or add to path.

[**IGV**](http://software.broadinstitute.org/software/igv/). this is a visualizer of bam files.

[**bedtools**](https://bedtools.readthedocs.io/en/latest/). Very powerful tool for manipulating and filtering bed files. You will need it at some point during your career.

**NOTE**: To add a given binary located in folder $DIRBIN to your path

	export PATH=$DIRBIN:$PATH
