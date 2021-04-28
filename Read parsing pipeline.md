[[Projects]]

Date:: [[October 5th, 2020]]

Due date:: 

Completed date:: 

Outcome:: Reusable modules that can be built into different pipelines

Success criteria:: 

Waiting for:: 

Status:: #started 

Priority:: #high

Tags:: #project #programming #gihub #pipeline 

Related projects:: [[ME/CFS]], [[Pipeline for UHT-ATAC]], [[Pipeline for sci-ATAC]]

Notes:: 

Next actions:: 

{{[[DONE]]}} Update notes from [[January 5th, 2021]] meeting #gihub #pipeline 

{{[[TODO]]}} Make notes on [[Jen Grenier]] pipeline document #gihub #pipeline #high 

{{[[TODO]]}} Run through [[ArchR]] tutorial #gihub #pipeline #high 

Did they write one for barnyard experiments? - No

{{[[DONE]]}} Get back to [[Jen Grenier]] on what mapping software do they recommend to produce tables for processing by [[ArchR]]? #gihub #pipeline 

{{[[TODO]]}} Look into [[UMI_tools]] #gihub #pipeline 

{{[[TODO]]}} Look into passing arguments to [[Shellscript]] using opt #gihub #pipeline 

{{[[TODO]]}} Make notes on [[cutadapt]] #gihub #pipeline #high 

{{[[TODO]]}} Look into [[MultiQC]] #gihub #pipeline 

{{[[TODO]]}} Look into [[STAR]] #gihub #pipeline 

{{[[TODO]]}} Look into [[ZUMI]] - we're not currently using this, but might be in the future #gihub #pipeline 



Producing barcodes.fa file from spreadsheet

Producing barcodes.fa file from spreadsheet

To run (from sandbox/GitHub/GG02/Read_parsing_pipeline/scripts directory):

`python3 produce_barcode_file.py ../../../../../read_parsing_pipeline/GIH-sciATACseq-2level-sampleDemux-BoExample.xlsx -o i7tagBC_i5tagBC_barcodes.fa`



[[March 30th, 2021]] [[Read parsing pipeline]]

[[March 30th, 2021]] [[Read parsing pipeline]]

Make the output file/logfile names more coherent/consistent

{{[[TODO]]}} What about demux.log file prefix? (line 175)

Working on using [[UMI_tools]] instead of [[Picard]]



[[March 29th, 2021]] [[Read parsing pipeline]]

[[March 29th, 2021]] [[Read parsing pipeline]]

{{[[TODO]]}} Question for [[Faraz Ahmed]]: passing barnyard genome to [[HOMER]] for the makeTagDirectory step - why did it work without this?

{{[[TODO]]}} multiQC - do we need the older version, or can we run the one BioHPC has installed?

{{[[TODO]]}} Look into SNAP-ATAC

{{[[TODO]]}} Need to add fastQC back into the pipeline

{{[[TODO]]}} Need to use UMI_tools instead of Picard



[[March 26th, 2021]] [[Read parsing pipeline]]

[[March 26th, 2021]] [[Read parsing pipeline]]

Building 100,000 record test sets:

2,041,600,588 rows in dataset - need to select records from around the middle

zcat CKDL210001568-1a_HFL3TCCX2_S1_L008_R1_001.fastq.gz | head -n 1030000000 | tail -n 400000 | gzip > CKDL210001568-1a_HFL3TCCX2_S1_L008_medium_R1_001.fastq.gz

zcat CKDL210001568-1a_HFL3TCCX2_S1_L008_R2_001.fastq.gz | head -n 1030000000 | tail -n 400000 | gzip > CKDL210001568-1a_HFL3TCCX2_S1_L008_medium_R2_001.fastq.gz

zcat CKDL210001568-1a_HFL3TCCX2_S1_L008_I1_001.fastq.gz | head -n 1030000000 | tail -n 400000 | gzip > CKDL210001568-1a_HFL3TCCX2_S1_L008_medium_I1_001.fastq.gz

zcat CKDL210001568-1a_HFL3TCCX2_S1_L008_I2_001.fastq.gz | head -n 1030000000 | tail -n 400000 | gzip > CKDL210001568-1a_HFL3TCCX2_S1_L008_medium_I2_001.fastq.gz

Running and testing horizontal merge in [[Pipeline for sci-ATAC]]



[[March 25th, 2021]] [[Read parsing pipeline]]

[[March 25th, 2021]] [[Read parsing pipeline]]

Running and testing horizontal merge in [[Pipeline for sci-ATAC]]



[[March 24th, 2021]] [[Read parsing pipeline]]

[[March 24th, 2021]] [[Read parsing pipeline]]

Producing barcodes.fa file from spreadsheet

To run (from sandbox/GitHub/GG02/Read_parsing_pipeline/scripts directory):

`python3 produce_barcode_file.py ../../../../../read_parsing_pipeline/GIH-sciATACseq-2level-sampleDemux-BoExample.xlsx -o i7tagBC_i5tagBC_barcodes.fa`



[[March 23rd, 2021]] [[Read parsing pipeline]]

[[March 23rd, 2021]] [[Read parsing pipeline]]

Modifying horizontal merge script per meeting with [[Jen Grenier]] on [[March 22nd, 2021]]



[[March 22nd, 2021]] [[Read parsing pipeline]]

[[March 22nd, 2021]] [[Read parsing pipeline]]

Things to talk to [[Jen Grenier]] about:

{{[[DONE]]}} What do we do with the quality scores on line 4 of the fastq record - do we pad with "I" for the extra barcode lengths?

Write the same $ variables as for line 2 - this will use the actual quality scores for the extra barcodes

{{[[DONE]]}} Do we want to set up /workdir/software for our own installations (e.g. [[HOMER]] and [[featureCounts]])

Yes

{{[[DONE]]}} Faraz running old version of multiQC - need to talk to him when he gets back

{{[[DONE]]}} Can't get atacQC (last step to run) - need permissions to a few things in [[Faraz Ahmed]] home directory

Setting up barnyard genome in [[HOMER]] configuration

Installing [[featureCounts]] in my workdir

i7tagBC is right before the ME sequence - 10 nt before the $3 position - 4 - 9 nt is the UMI, 10nt i7tabBC, and then the 19nt ME sequence

Need to write code to parse the whitelist specific to a project, grab the sample names, and write these out to an .fa file that [[cutadapt]] can use as input

Modifying horizontal merge script per today's meeting



[[March 19th, 2021]] [[Read parsing pipeline]]

[[March 19th, 2021]] [[Read parsing pipeline]]

Installed [[featureCounts]] in /workdir/subread/ directory

[[HOMER]] needs to be passed the barnyard genome - need to check with [[Faraz Ahmed]] on this

To configure a custom genome:

Website: http://homer.ucsd.edu/homer/introduction/update.html



[[March 18th, 2021]] [[Read parsing pipeline]]

[[March 18th, 2021]] [[Read parsing pipeline]]

Continuing to modify horizontal merge script



[[March 17th, 2021]] [[Read parsing pipeline]]

[[March 17th, 2021]] [[Read parsing pipeline]]

Modifying horizontal merge script to add duplicate barcodes to 5' end of the R2 read:

https://stackoverflow.com/questions/15969861/pass-parameter-to-an-awk-script-file

https://unix.stackexchange.com/questions/475008/command-line-argument-in-awk

Installing [[HOMER]] v4.11.1 in /workdir/prm88/homer directory

Note: [[Faraz Ahmed]] was using v4.10.4





[[March 16th, 2021]] [[Read parsing pipeline]]

[[March 16th, 2021]] [[Read parsing pipeline]]

Building [[Pipeline for sci-ATAC]]

Getting a better idea of how i5 and i7 barcodes uniquely identify a sample

Looking into using [[UMI_tools]] to do the markDups step, instead of Picard



[[March 15th, 2021]] [[Read parsing pipeline]]

[[March 15th, 2021]] [[Read parsing pipeline]]

Building [[Pipeline for sci-ATAC]] - adding link to human / mouse combined genome ref

BWA manual: http://bio-bwa.sourceforge.net/bwa.shtml

Meeting notes with [[Jen Grenier]], [[Adrian McNairn]], and [[Elizabeth Fogarty]]

![](https://firebasestorage.googleapis.com/v0/b/firescript-577a2.appspot.com/o/imgs%2Fapp%2FWorkingmemory%2FW3pnVb9LL9.png?alt=media&token=20cf4ef0-be3d-44ef-a943-c8eccf4a4ca1)

Demultiplexing is usually at the sample level, using the index reads that come from the PCR

The reason we need to demultiplex on the tagmentation component of the barcode is to maintain unique barcodes (normally people use PCR)

However, we have a limited number of PCR barcodes, and we need to use all of our PCR barcodes in every sample.

So we have more tagmentation barcodes than we have PCR barcodes, so we use different sets of tagmentation barcodes for different samples and use those to demultiplex a shared lane

However, when demultiplexing on the i5 tagmentation barcode, we also need to keep that sequence in the fastq file. Since cutadapt will remove the i5 barcode, we need it to be duplicated so that one copy remains after running cutadapt

Questions for Bo:

![](https://firebasestorage.googleapis.com/v0/b/firescript-577a2.appspot.com/o/imgs%2Fapp%2FWorkingmemory%2FZhH6U45YJw.png?alt=media&token=56f8b7b2-f07a-4efe-9f13-74fd6244cbb1)

Give us a barcode library in the following format: xxx

What are the barcode sequences and names?

As the programmer, what format would I like to see?



[[March 12th, 2021]] [[Read parsing pipeline]]

[[March 12th, 2021]] [[Read parsing pipeline]]

Building [[Pipeline for sci-ATAC]]

Downloaded data from [[Paul Soloway]] and [[Blaine Harlan]]



[[March 11th, 2021]] [[Read parsing pipeline]]

[[March 11th, 2021]] [[Read parsing pipeline]]

Modifying / testing [[Pipeline for UHT-ATAC]]

Building [[Pipeline for sci-ATAC]]



[[March 10th, 2021]] [[Read parsing pipeline]] Modifying [[Pipeline for UHT-ATAC]]

[[March 10th, 2021]] [[Read parsing pipeline]] Modifying [[Pipeline for UHT-ATAC]]



[[March 9th, 2021]] [[Read parsing pipeline]] Modifying [[Pipeline for UHT-ATAC]]

[[March 9th, 2021]] [[Read parsing pipeline]] Modifying [[Pipeline for UHT-ATAC]]



[[March 8th, 2021]] [[ME/CFS]] weekly meeting

[[Read parsing pipeline]]

Intermediate demultiplexing needed, using the i5 barcodes

These are sample specific barcodes

Might have to do this after the horizontal merge step, because of I1 and I2 (current demultiplexing methods only work for two files, R1 and R2)

For [[Adrian McNairn]] last experiment, why are there only 600 barcodes? Should be 2 - 5,000

Barcode whitelists for our experiments are stored on the server - each experiment uses a subset of these

![](https://firebasestorage.googleapis.com/v0/b/firescript-577a2.appspot.com/o/imgs%2Fapp%2FWorkingmemory%2FiQts9bCAoC.png?alt=media&token=94f2fcd0-c2a0-4277-90d7-f4a61362235a)



[[March 8th, 2021]] [[Read parsing pipeline]] Modifying UHT-ATAC pipeline

[[March 8th, 2021]] [[Read parsing pipeline]] Modifying UHT-ATAC pipeline



[[March 5th, 2021]] [[Read parsing pipeline]] Modifying UHT-ATAC pipeline

[[March 5th, 2021]] [[Read parsing pipeline]] Modifying UHT-ATAC pipeline



[[March 4th, 2021]] [[Read parsing pipeline]] Making notes on existing UHT-ATAC pipeline

[[March 4th, 2021]] [[Read parsing pipeline]] Making notes on existing UHT-ATAC pipeline

[[Pipeline for UHT-ATAC]]



[[March 3rd, 2021]] [[Read parsing pipeline]] Making notes on yesterday's meeting and [[UMI_tools]]

[[March 3rd, 2021]] [[Read parsing pipeline]] Making notes on yesterday's meeting and [[UMI_tools]]



[[March 2nd, 2021]] [[Read parsing pipeline]] meeting with [[Jen Grenier]] and [[Faraz Ahmed]]

[[March 2nd, 2021]] [[Read parsing pipeline]] meeting with [[Jen Grenier]] and [[Faraz Ahmed]]

[[Jen Grenier]] had asked (yesterday) about the organization of files on the Danko cluster, and was that a structure we should consider for our server:

Divided into data and projects

Data is non-project specific datasets such as genomes (divided by species, version, cell type), and other annotations (PRO-seq tracks, GRO-seq tracks, histone marks, ChIP-seq tracks, poly-A tracks, etc.)

Projects is divided first by net ID and then by individual project - contains datasets such as fastq files

Links from [[Faraz Ahmed]] following the meeting:

UMI tools dedup: https://umi-tools.readthedocs.io/en/latest/reference/dedup.html

Blog giving reason for parameter values used when ATAC-seq peak calling with [[MACS]]: https://www.biostars.org/p/209592/

More on peak calling for ATAC-seq: https://bioinformatics-core-shared-training.github.io/cruk-autumn-school-2017/ChIP/Materials/Lectures/Lecture5_Peak%20Calling_SS.pdf

Video of call: C:\Users\prm88\Documents\Box\Zoom\2021-03-02 15.46.34 TREx_Jen Grenier_ Office Hours 658238586

Started with GI in-house analysis plans.pptx

Showing UHT-ATAC-seq v0.2

![](https://firebasestorage.googleapis.com/v0/b/firescript-577a2.appspot.com/o/imgs%2Fapp%2FWorkingmemory%2FZngZAdk5hk.png?alt=media&token=f4b7f76a-8675-4bbe-9252-ebc626efda3c)

This is high throughput ATAC, so it uses sample bar codes

Page 2 shows sci-ATAC-seq (2 level)

![](https://firebasestorage.googleapis.com/v0/b/firescript-577a2.appspot.com/o/imgs%2Fapp%2FWorkingmemory%2F3MePLmmTAn.png?alt=media&token=b26ee6b6-7855-42c1-9d32-97fdad9f4886)

The first part is pretty non-standard and requires some custom shell script: `horizontalMerge-padUMI.sh`

This results in a pair of fastq files, where the fastq headers have the cell bar code and the UMI.

Following this we have the standard trimming of adapters and filtering for quality and short reads

Should be ready for mapping after this step

All code currently in [[Jen Grenier]] directory on GG02

GenomicsInnovation/sciATAC-2level/GIA0130-AdrianBarnyard/fastqs/ 

Raw fastq files

GenomicsInnovation/sciATAC-2level/GIA0130-AdrianBarnyard/fastqs/parsed-BCs

Parsed bar code files

GenomicsInnovation/sciATAC-2level/GIA0130-AdrianBarnyard/fastqs/parsed-BCs/BCwhitelist_moreRuns

Might need one more trimming step after this

The whitelist involved at this step catches single nucleotide errors in the barcode



For the bulk ATAC-seq project from last August: (12:00)

GenomicsInnovation/UHT-ATAC

A lot of cross talk

In this case the reads have been mapped - in the bwa directory

Uses a shell script from [[Faraz Ahmed]] to run bwa (bwa-mem.sh) (12:39)

![](https://firebasestorage.googleapis.com/v0/b/firescript-577a2.appspot.com/o/imgs%2Fapp%2FWorkingmemory%2FMjaokfwJN-.png?alt=media&token=c34e03c0-ed20-4f3a-96ba-a7d9e2196ba2)

Takes advantage of 10x's hg39 / mm10 combined reference

Jen's sorted files are quite a bit smaller than the unsorted - might have to resort in case they are corrupted (Jen had problems using a pipe in the sort and had to do it in two steps)

Per Faraz's recommendation, don't map and sort at the same time

Separate scripts for sorting: samtool-sort.sh and samtools-sort-stats.sh



Switched over to Faraz's pipeline

If you put the processed files in a folder called trim, the entire pipeline should run

The trimPE and trimHiSeqPE functions are essentially the same, but trimPE handles data from NextSeq and trimHiSeqPE handles data from NovaSeq

In both of these, length should be a variable, and set the default to 50

For the alignPE function, pretty much the same as the code Jen was showing above

The next step, sort, sorts the BAM file and gathers stats for the QC

rmMT removes mitochondrial reads

Next step: mark and remove duplicates

This uses picard tools, but we should try it using UMI tools - probably more applicable to the files we have

Then we generate more stats and run multiqc

The multiQC will use all stats files, so the ones generated at this step and the ones generated above

Next we use HOMER to generate the bedgraph file that we use for visualization

Next we do peak calling using MACS2

In the modular system we're designing, may want to allow users to select which annotation / peak calling method they use

Blog giving reason for parameter values used when ATAC-seq peak calling with [[MACS]]: https://www.biostars.org/p/209592/

We call peaks on each individual BAM file and then merge peaks for all peak calls combined

Merged narrowPeak file is then used to create saf file

{{[[TODO]]}} What is a saf file?

Use saf file with feature counts to assign each fragment to a given merged peaks site

{{[[TODO]]}} What is a frip file?

Then we produce the bedgraph files

Initially we need to answer questions like what are the relative amounts of mouse vs. human

Do the frip calls get better or worse when you change the temperature

{{[[TODO]]}} What was the name of the tool that [[William Lai]] mentioned - might be better at handling ATAC-seq data





[[January 29th, 2021]] Meeting with [[Jen Grenier]] and [[Blaine Harlan]]

[[January 29th, 2021]] Meeting with [[Jen Grenier]] and [[Blaine Harlan]]

Related projects:: [[Read parsing pipeline]]

{{[[DONE]]}} Write up notes on meeting with [[Blaine Harlan]] #gihub 

Used one row / lane

~200ng / sample

multiqc after sequencing

demultiplex - trim for adapter and low quality bases

gave Blaine barcodes

Problem: running UMI dedup gives really high high rate of duplication

[[DESeq]] should handle the lower detection rates

Remove low gene counts prior to running [[DESeq]]

Make histogram of normalized counts / gene and use that to decide where cutoff should be. Should give a normal curve = cutoff before it gets 'noisy.'

Top 500 variable genes drive the PCA plots

{{[[TODO]]}} How would I figure this out in [[R]] / [[Seurat]] / [[ArchR]]? #gihub  



[[January 12th, 2021]] [[Read parsing pipeline]] meeting with [[Jen Grenier]] and [[Faraz Ahmed]]

[[January 12th, 2021]] [[Read parsing pipeline]] meeting with [[Jen Grenier]] and [[Faraz Ahmed]]

Video of call: C:\Users\prm88\Documents\zoom\2021-01-12 15.18.06 TREx_Jen Grenier_ Office Hours 658238586

UHT-3' RNA-seq

![](https://firebasestorage.googleapis.com/v0/b/firescript-577a2.appspot.com/o/imgs%2Fapp%2FWorkingmemory%2FDBh8w0rcdq.png?alt=media&token=96c97427-2e83-46a9-98f9-eee9219b75b2)

Architecture of library:

The arrows show were the primers are that we get off of the sequencer (up to 4)

Primer followed by the sequence of mRNA

Read 2 is on the opposite strand going in the other direction, so it's the paired end in the insert read

The first thing read 2 generates is a UMI sequence, then an inline bar code, then a poly-T stretch, and then it gets into the insert from the other end

It also generates an index read that could be at the beginning as a PCR bar code.

These index reads are used to demultiplex by sample

So the PCR bar code are bar coding pools of samples

If the run is not demultiplexed you get 4 FASTQ files (read1, read2, Index1, and Index2) - each file contains all clusters on the flowcell

Critical that we keep the rows in each file in sync with each other

---

Step 1: Move UMI to R1/R2 headers as early as possible (when everything is still in sync)

The first 8 characters in read 1 and read 2

Read 2 contains the UMI

Step 1B: Horizontal merge I1 and I2 reads into R2

These next steps are about custom demultiplexing, using the bar code information

There may already be some demultiplexing that happened on the PCR bar codes

Note: it's called an inline bar code because it's in the insert read, not in the index read

Used paste command to create longer inline bar codes at he 5' end of R2

Step 1B is optional for the example [[Jen Grenier]] was talking about, but will be required for combinatorial single-cell.

At this step, we've collapsed down to two files (14:00)

At this point we can use [[cutadapt]] to do the demultiplexing (14:50)

The goal here (for bulk UHT 3' RNA-seq) is to generate a pair of FASTQ files for each sample, with UMI's in the headed to take advantage of the combinatorial bar code (15:20)

After the paste, read 1 will be longer than read 2 (we moved the UMI out of read 2, so it's down to 142 characters)

Haven't yet tested whether it's ok to have Read 1 and Read 2 different lengths (22:20)

Another variable is the length of the index read - the standard is 8, but custom pipelines can change this

[[cutadapt]] allows you to provide a human readable name to go with each of the bar codes in your whitelist

![](https://firebasestorage.googleapis.com/v0/b/firescript-577a2.appspot.com/o/imgs%2Fapp%2FWorkingmemory%2F9lUp_HeAj_.png?alt=media&token=2a7f6dc0-a61a-403a-879d-e4d7a1a1076d)

Possible variables used in script:

UMI = 8 nt (removed)

Output = 142 nt (assumes R2 starts at 150 nt)

[[cutadapt]]

Focus on trimming paired-end reads

Keeping reads in sync

Demultiplexing

Example read qc in Slack

Prefer BWA over STAR for alignment of ATAC-seq libraries (1:15:00)

SAM tools now more strict some problems working with BWA - ask [[Faraz Ahmed]] about this

BWA only writes SAM files - by using regroup you can add a lot of stuff to the headerr

Use -b to produce a BAM file





[[January 5th, 2021]] [[Read parsing pipeline]] Meeting with [[Jen Grenier]]

[[January 5th, 2021]] [[Read parsing pipeline]] Meeting with [[Jen Grenier]]

Location of [[Jen Grenier]] notes:

/workdir/jkg47/GenomicsInnovation/UHT-ATAC/Novogene-G201123J-noDemux/fastqs

{{[[TODO]]}} Seurat vignette: https://satijalab.org/seurat/v3.1/de_vignette.html #gihub #pipeline 

Video of call: C:\Users\prm88\Documents\Zoom\2021-01-05 15.06.46 TREx_Jen Grenier_ Office Hours 658238586



[[UHT-ATACseq]] being done by [[Adrian McNairn]]

Analysis pipeline:

[[sciATACseq]]

Combinatorial barcode: i7 barcode + i5 tag + i7 tag

100 samples; 96 wells

{{[[TODO]]}} Relevant paper: Cusanovich (2015 or 2017) #gihub #pipeline 

i7 barcode embedded in read2 



Test dataset: GG02/workdir/jkg47/genomicsinnovationhub/raw_data

DNA barcodes package - design of barcodes



{{[[TODO]]}} Look into [[Case]] #gihub #pipeline 

{{[[TODO]]}} Look into [[Esac]] #gihub #pipeline 

{{[[TODO]]}} Look into [[Snakemake]] #gihub #pipeline 

{{[[TODO]]}} Look into [[Nexflow]], [[NextflowWorkbench: Reproducible and Reusable Workflows forBeginners and Experts]] #gihub #pipeline 

{{[[TODO]]}} Look into [[STARsolo]] - mapping, demultiplexing and quantification for single cell RNA-seq #gihub #pipeline 

{{[[TODO]]}} Look into [[split pool]] processing #gihub #pipeline 

{{[[TODO]]}} Look into [[Signac]] #gihub #pipeline 

{{[[TODO]]}} Look into [[Simpson Index]] #gihub #pipeline 



Email from [[Jen Grenier]]

I’ve attached two slides on my first attempts at preprocessing fastq files (parsing UMI and BCs) for the two UHT methods being developed in the GIH lab.

Note that the steps for [[UHT-3’Seq]] could be reordered to match [[UHT-ATACseq]], and allow for use of the long R2 read into the insert. In the test data, the R2 quality was really low so I didn’t try to use it, I just truncated to use only the UMI+RTBC. But to preserve the PE insert read from R2, it will work better to move the [[UMI seq]] first and then demux on RTBC that becomes 5’anchored in R2 (as for [[UHT-ATACseq]]).

We have a test dataset for each method on GG02:

/workdir/jkg47/GenomicsInnovation/UHT-3pSeq/Novogene-G200831J/

/workdir/jkg47/GenomicsInnovation/UHT-ATAC/Novogene-G200831A-rptdemux/

Feel free to copy any files and try to read through my notes on the server!

And of course ask any questions that come up.

For example, you may want to know what BC combos are expected – I will work on a ‘sample sheet’ that describes what we did in the lab for these test experiments.

Slides for [[UHT-ATACseq]] and [[UHT-3’Seq]]: D:\Box Sync\genome_innovation_hub\read_parsing_pipeline



