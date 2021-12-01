# NGS Crash Course
A few hour introductory course to NGS analysis, mainly for variant calling.

This is a few hour lecture that I usually teach at [Universitat Autonoma of Barcelona](www.uab.es) [MSc in Bioinformatics](https://mscbioinformatics.uab.cat/). It contains first basic steps for someone who has not had exposure to NGS analyses, yet is familiar with linux commands. For simplicity, I base teaching in one of the smallest microorganism ([Mycoplasma genitalium](https://www.ncbi.nlm.nih.gov/genome/?term=Mycoplasma%20genitalium)) and a few simulated reads using [ART](https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm). 

There exists dozens of procedures specific to many of the applications of NGS data.

**These notes are intended for variant calling with genome data.**

## WARNINGGGGG 1: This is a rapidly evolving field. 

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

### WARNINGGGGGS 2
- **NGS data are massive yet VERY NOISY.** 
- **You will be puzzled about how results can change depending on bioinformatics options utilized.**
- **Do NOT outsource this step. DO IT YOURSELF or someone you know well.**
- **Caution and quality control are a must in every step, especially when analyzing several samples simultaneously.**
- **Missing data are unavoidable.**

## QUALITY CONTROL IS A MUST IN EVERY NGS ANALYSIS STEP 
### VISUALIZATION TOOLS CAN HELP.

### What techniques are there?
Classical Sanger sequencing was a 'First' generation technology. Second or next generation sequencing was initially [454 technique](https://en.wikipedia.org/wiki/454_Life_Sciences) quickly supersed by [Illumina's technologies](https://en.wikipedia.org/wiki/Illumina_dye_sequencing). More recently, the third generation sequencing technologies deliver much longer reads in lower quantity and of lower quality. Representatives are [PacBio](https://www.sciencedirect.com/science/article/pii/S1672022915001345) and [Oxford nanopore](https://en.wikipedia.org/wiki/Nanopore_sequencing).

Each of these technologies may require slightly different software / options but bioinformatics pipelines are essentially the same.

### Which are the main steps?
> 1. Sequence quality check is warranted using [fastq](https://en.wikipedia.org/wiki/FASTQ_format) or similar
> 2. Reference genome indexing.
> 3. Sequence alignment against a reference genome.  
> 4. Alignment file needs to be sorted and can be polished by removing PCR duplicates, quality can be recalibrated with a known set of variants.
> 5. Variant calling and quality filtering.

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

### Genome indexing
- It consists of generating a rapidly accessible reference genome. 
- It needs to be done only one.
- Currently,  [Burrows-Wheeler algorithm](https://en.wikipedia.org/wiki/Burrows%E2%80%93Wheeler_transform) is the most widely used 
- BWA software is the most used software ([Li and Durbin, 2019, Bioinformatics](https://academic.oup.com/bioinformatics/article/25/14/1754/225615)).
- It can be slow but needs to be done only once.

The instruction in bwa is:

	# assembly.fa is a fasta file with reference genome
	# type bwa index for options available
	bwa index -a bwtsw assembly.fa

### Reads alignment
* It consists of determining the most likely origin of a set of short read sequence within a larger reference genome.
* BLAST / CLUSTALW are classical tools, but you will never finish …
* New challenges: number of sequences (speed) and close similarity (large impact of sequence errors).
* Main Softwares: 
	- [BWA](https://sourceforge.net/projects/bio-bwa/files/) (genome data)
	- [Hisat](http://daehwankimlab.github.io/hisat2/hisat-3n/) (RNAseq)	- 
* NGS read alignment is not an error free process nor is 100% error free. Errors can occur because wrong or incomplete genome references, incomplete search, indels, gaps, segmental duplications,...

### Variant calling
* It is the main goal of many studies.
* Fraught with dangers and subtleties, among them:
	- Base and mapping qualities
	- Low or too high coverage
	- Multiple alignments
	- ...
* Structural variants are much more difficult to identify than SNPs because reads are short and reference genomes usually contain numerous errors.
* Sex chromosomes require specific algorithms
* Pools also require specific software (e.g., [Raineri et al. 2012](https://doi.org/10.1186/1471-2105-13-239)).
* The standard format to contain SNPs is the [vcf format](https://github.com/samtools/hts-specs)

There are essential differences between how SNPs are identified with Sanger and NGS sequencing. 

![](https://github.com/miguelperezenciso/NGScrashcourse/blob/master/sanger.png)

In Sanger sequencing (Figure above), a polymorphic site is identified by similar intensities of each nucelotide fluorescence. It is usually quite accurate procedure.

![](https://github.com/miguelperezenciso/NGScrashcourse/blob/master/igv.png)

With NGS data, variant calling is done by comparing reads against a reference a genome, a putative SNP is identified by counting the number of reads that have the reference vs those that have an alternative allele, base qualities is also considered. Since read depth is stochastic and varies largely along the genome, there will always be regions of the genome sequenced with low depth. In these regions, it is a matter of how much risk you are ready to take on whether you take a SNP as real or not. You will encounter this dilemma even at high read depths. In conclusion, SNP calling is a probabilistic decision. **For these and other reasons, raw variant calls need to be filtered.** 

One parameter you need to filter for is read depth. Variants from very low depth (say < 5 - 10 depending on your average read depth) can be removed. Variants from very high depth are equally unreliable, as they can pertain to repeated regions. I recommend to filter out variants from regions with > twice average depth.

**Some advices:**
- Using several individuals simultaneously for SNP calling improves reliability, specially for middle frequency alleles.
- Distrust indels, especially long or complex ones.
- It is absolutely necessary to filter SNPs with bcftools filter or similar.
- Inspect with IGV if necessary.

### WARNINGGGG: Variant filtering is essential to improve reliability

### Visualization
By its own nature, NGS data analyses need to be done with automatic tools and is impossible to visualize each result, eg, a single mmamalian genome can contain millions of SNPs. However, visualizings some specific regions can be recommended. The most popular tool is IGV. It requires an indexed bam file and a reference genome.

### Exercise
![](https://github.com/miguelperezenciso/NGScrashcourse/blob/master/Xchr.png)

	This figure represents the average depth of an X chromosome from a mammalian individual. Is the individual sequenced a male or a female? Why?

### SNP annotation
By annotation we mean identifying general genome features associated with a SNP or a genome region. For instance, SNP annotation means identifying whether the SNP resides in a coding, regulatory, intergenic region, if coding whether it is synonymous or not, if potentially deleterious or not, and so on. The easiest to use tool is Ensembl Variant Effect Predictor ([VEP](https://www.ensembl.org/info/docs/tools/vep/index.html)). 

### Experimental design
Obtaining sequence is far more expensive than genotyping using a commercial array. Sequencing provides a list of variants free of ascertainment bias, as opposed to array genotyping. Choice of technology then needs to be driven by the target of the experiment. Even if sequencing allows to uncover all variants, it is not free from limitations: (i) expensive, (ii) cumbersome bioinformatics analysis, (iii) variant calling is an error prone process itself. For a given budget, how would you allocate resources, would you sequence more individuals at sahllow depth or few at high depth?

We explored this issue in [Nevado et al. (2014a)]( https://doi.org/10.1111/mec.12693) and in [Nevado et al. 2014b](https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.12286). We find that relying on genotype calls provides biased estimates of population genetic statistics at low to moderate read depth (2–8×). Direct estimation without calling genotypes returns the most accurate estimates of variability and of most SFS tests investigated, including at low read depth (2–4×). Studies without species-specific reference genome should aim for low read depth and avoid genotype calling whenever individual genotypes are not essential. Otherwise, aiming for moderate to high depth at the expense of number of individuals is recommended.

![](https://github.com/miguelperezenciso/NGScrashcourse/blob/master/pipeliner.png)

The above figure from nevado et al. (2014b) shows several strategies with similar costs, barred from library contruction costs, that range from 50 individuals at 4x depth to 4 individuals at 50x depth. Upper panel shows the percentage of original heterozygous SNP genotypes that are correctly identified and bottom panel shows the number segregating sites identified with each experimental design. Note a depth of 10 - 20 x can be an optimum, increasing beyond that does not increase recovery but decreases the number of sites found because fewer individuals are sequenced.

### RNAseq
Sequencing RNA vs expression arrays offers some important advantages:
* Discovery of new transcripts, isoforms, and non coding RNAs.
* Quantification of (differential) expression. 
* Allele specific expression, 
* Improve anotation.

Best if annotation file available.

Paired end needed if you wish to discover isoforms, SE is enough if only expression quantification is required.

Mapping is normally done allowing ambiguity to study expression in paralogs.

Popular pipeline: [Hisat2](http://daehwankimlab.github.io/hisat2/).


### Metagenomics

# TAKE HOME MESSAGES
- NGS data are massive yet noisy. Quality check and filtering is a must in every step. Visualizing bam tools can be useful.
- Experimental design is relevant. Is NGS really what you need? 
- Numerous pipelines and software are available, specific to each application (Variant detection, RNAseq, metagenomics, …).
- Analyses are full of subleties, do the analysis yourself.
- Allow enough computer power, cloud services are good options.

==========
# PIPELINE
==========

### Folder scheme
I suggest to have separate folders to organize the different analysis steps, but this is very personal. In the exercise, we have a folder **[assembly](https://github.com/miguelperezenciso/NGScrashcourse/tree/master/assembly)** with the assembly and all required indices, a folder **[reads](https://github.com/miguelperezenciso/NGScrashcourse/tree/master/reads)** with all sequence data, bam files are stored in **bamfiles** directoy and vcf files in **varfiles**. Optionally, a **[bin]((https://github.com/miguelperezenciso/NGScrashcourse/tree/master/bin)** folder contains executables. Alternatively, you can have links to these tools in your main path.

	# current folder
	DWD=$(pwd)
	# this folder should contain the assembly
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

### Main tools needed
Dozens of tools have been developed for NGS data, and each kind of data (eg, genomic, RNAseq, microbiome...) requires its own specific tools. Here I focus on SNP calling. 

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

### EXERCISE 2: 
Install these softwares in the **bin** directory. Note: IGV can be run via [web](https://igv.org/app/). We will not use bedtools here.

### Definitions
Here I assign values to some variables such that the pipleine can be easily changed. Remmember in shell $NAME refers to the value of variable NAME.

	# define variables with program names or install in your path
	fastqc=$DIRBIN/FastQC/fastqc
	bwa=$DIRBIN/bwa
	samtools=$DIRBIN/samtools
	bcftools=$DIRBIN/bcftools
	gatk=$DIRBIN/gatk-4.2.3.0/gatk  # adjust version as needed
	vcftools=$DIRBIN/vcftools

	# sample to be analyzed
	OUT=Individual1

	# reads from sample $OUT, they are in $DIRDATA
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
   	$gatk CreateSequenceDictionary R=$DIRASSEMBLY/$ASSEMBLY.fna O=$DIRASSEMBLY/$ASSEMBLY.dict
   
	# back to working directory
	cd $DWD

### Step 1: Check reads quality

	cd $DIRDATA
	zcat $DIRDATA/$READS_PE1 | $fastqc stdin
	zcat $DIRDATA/$READS_PE1 | $fastqc stdin
	# move back to wkng folder
	cd $DWD

You should see something like 

![](https://github.com/miguelperezenciso/NGScrashcourse/blob/master/fastqc.png)

### Step 2 BWA alignment and refinment
This is the most expensive part, both in CPU and memory usage

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

   	# remove PCR duplicates
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

   	# back home
	cd $DWD

### Step 3a. samtools single sample SNP CALLING with gVCF blocks 
This is the trickiest step, the one with most options and highly dependent on depth and sequence quality. Filtering SNPs is a must. Distrust indels more than SNPs. Beware of options with snp calling using several samples together. In this step we use samtools / bcftools pipeline.

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
	
### Step 3b. GATK Multiple sample SNP CALLING with gVCF blocks
Based in https://gatk.broadinstitute.org/hc/en-us/articles/360035890411?id=3893

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

### Step 4. Visualize some SNPs with IGV (http://software.broadinstitute.org/software/igv/) 

	# Start IGV (you need java8)
	# You need to load the Micoplasma genome in $DIRASSEMBLY
	java -Xmx1000m -jar $igv
	
You can also run via web app: https://igv.org/app/

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
