# NGScrashcourse
A few hour introductory course to NGS analysis

This is a few hour lecture approx that I teach at Universitat Autonoma of Barcelona MSc in Bioinformatics. It contains first basic steps for someone who has not had exposure to NGS analyses, yet is familiar with linux commands. For simplicity, I base teaching in one of the smallest microorganism (Mycoplasma genitalium) and a few simulated reads using [ART](https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm).

### Exploring...
1- Browse GTAK best practices (https://software.broadinstitute.org/gatk/best-practices/):In particular, you can take a look at  https://software.broadinstitute.org/gatk/best-practices/workflow?id=11165; and https://software.broadinstitute.org/gatk/best-practices/workflow?id=11145

    
2- Browse the SRA archive (https://www.ncbi.nlm.nih.gov/sra/docs), which holds NGS data of numerous types and species. Project ERP004545 corresponds to M genitalium.

## Exercises
1 Choose a different strain from experiment ERP004545, you need to install aspera

2 Run fastqc, Perform alignment as described

3 Perform multiple snp calling with about 5-10 samples together

4 Using vcftools, do a plot of allele frequency and missingness across all SNPs
