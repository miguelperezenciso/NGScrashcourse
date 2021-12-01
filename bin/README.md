# INSTALLING TOOLS

Dozens of tools have been developed for NGS data, and each kind of data (eg, genomic, RNAseq, microbiome...) requires its own specific tools. Here I focus on SNP calling. 

![](https://github.com/miguelperezenciso/NGScrashcourse/blob/master/software.png)

You need to install:
[**FastQC**](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/): This is a java tool that provides graphics of sequence quality. 
Documentation can be found [here](https://raw.githubusercontent.com/s-andrews/FastQC/master/INSTALL.txt).
You can run FastQC in one of two modes, either as an interactive graphical application.

`
# running interactively in linux
./fastqc
# and navigate to load a compressed fastq file

# run as shell
zcat *fastq.gz | fastqc stdin
`

- [BWA](https://sourceforge.net/projects/bio-bwa/files/). This is one of the most widely used aligner for genomic data.

	`git clone https://github.com/lh3/bwa.git`
	`cd bwa; make` 

- samtools, bcftools and hstslib: http://www.htslib.org/download/. A series of tools for manipulating bam files and variant calling. Once extracting, for each of them type 'make' and move executables to bin or to add to path. 

- [GATK](https://gatk.broadinstitute.org/hc/en-us).

- [vcftools](https://sourceforge.net/projects/vcftools/). Filters and extract info from vcf files, format conversion with plink. Download, do 'make' and move executable to bin or to add to path.

- [igv](http://software.broadinstitute.org/software/igv/). this is a visualizer of bam files.

- [bedtools](https://bedtools.readthedocs.io/en/latest/). Check also bedtools, you will need it at some point during your career.
