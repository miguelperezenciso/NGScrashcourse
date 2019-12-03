# NGScrashcourse
A few hour introductory course to NGS analysis

This is a few hour lecture that I teach at Universitat Autonoma of Barcelona [MSc in Bioinformatics](https://mscbioinformatics.uab.cat/). It contains first basic steps for someone who has not had exposure to NGS analyses, yet is familiar with linux commands. For simplicity, I base teaching in one of the smallest microorganism ([Mycoplasma genitalium](https://www.ncbi.nlm.nih.gov/genome/?term=Mycoplasma%20genitalium)) and a few simulated reads using [ART](https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm). 

### Folder scheme
I suggest to have separate folders to organize the different analysis steps, but this is very personal. In the exercise, we have a folder **assembly** with the assembly and all required indices, a folder **reads** with all sequence data, bam files are stored in **bamfiles** directoy and vcf files in **varfiles**. Optionally, a **bin** folder contains executables.

### Tools needed
Hundreds of tools have been developed for NGS data, and each kind of data (eg, genomic, RNAseq, microbiome...) requires its own specific tools. It focuses on SNP calling. I may add new pipelines in the future for other type of analyses. 

Note: This is a samtools based pipeline, alternative pipelines are based on GATK. Hopefully, GATK pipelines will be added soon.

You need to install:
- FastQC: 

### Exploring...
 - Get acquainted with major sequenicng technologies: Illumina (https://en.wikipedia.org/wiki/Illumina_dye_sequencing), Oxford Nanopore (https://en.wikipedia.org/wiki/Nanopore_sequencing), PacBio(https://en.wikipedia.org/wiki/Single-molecule_real-time_sequencing)...

 - Browse [GTAK](https://software.broadinstitute.org/gatk) best practices (https://software.broadinstitute.org/gatk/best-practices/):In particular, you can take a look at  https://software.broadinstitute.org/gatk/best-practices/workflow?id=11165; and https://software.broadinstitute.org/gatk/best-practices/workflow?id=11145
    
 - Browse the SRA archive (https://www.ncbi.nlm.nih.gov/sra/docs), which holds NGS data of numerous types and species. Project ERP004545 corresponds to M genitalium.

## Exercises
 * Choose a different strain from experiment ERP004545, you need to install aspera

 * Run fastqc, Perform alignment as described

 * Perform multiple snp calling with about 5-10 samples together

4 Using vcftools, do a plot of allele frequency and missingness across all SNPs
