## Description 

**Seq2Fun** is an ultra-fast, all-in-one functional profiling tool for RNA-seq data analysis for organisms without reference genomes.
Visit Seq2Fun websit **www.seq2fun.ca** for more details.


## Key features of Seq2Fun
Ultra-fast: Seq2Fun is > 120 times faster (~ 2 million reads / minute) than the conventional RNA-seq workflow.

Extremely low memory cost: Seq2Fun consumes as little as 2.27 GB memory and can run on a standard PC with 8 threads and 16 GB memory.

Highly efficient: Seq2Fun can finish a typical RNA-seq dataset within several hours on a standard PC in stead of several days or even weeks by conventional RNA-seq workflow on a high-performance server.

Reference-free: Seq2Fun does not require the genome or transcriptome reference of the organism; it is also transcriptome de novo assembly-free.

Highly accurate: Seq2Fun generates KO abundance with R2 value as high as 0.93 comparing with the ground truth.

All-in-one: Seq2Fun directly takes raw RNA-seq reads as input and output gene abundance table without any intermediate file writing and loading, making I/O very efficient.

Multifunctional: Seq2Fun generates 6 levels of output files, including KO abundance table, hit pathway table, hit species table, reads KO table, a html report summarizing these tables and reads quality check, as well as output mapped clean reads for further analysis such as gene assembly.

Highly flexible: Seq2Fun supports RNA-seq analysis on particular genes or groups of organisms using customized database.

Easy to use: Seq2Fun requires minimal programing skills.

## Getting started
### Step 1. Install the package
Seq2Fun (version 1.0.0) is written in C/C++11 and can be installed on Linux or Mac OS X (with Xcode and Xcode Command Line Tools installed). 
We have tested Seq2Fun on Ubuntu (16.04 LTS and above) and macOS Catalina.

Click [here](). to download the source code. From within the folder containing the downloaded package, issue the following commands:

```
tar -xvzf seq2fun_v1.0.0.tar.gz
Or 
git clone https://github.com/xia-lab/Seq2Fun.git

cd seq2fun/src/
make clean
make
```
### Step 2. Run a small test 
There are four sub folders under seq2fun - src, bin, database and testdata. The bin folder contains the binary code we just complied. The testdata contains a small test data from the Case Study.
From within the **testdata folder**, issue the following commands:
```
../bin/seq2fun --sampletable sample.txt --tfmi birds_cdhit99_proteins.fmi --genemap birds_protein_ko_species_cdhit99.txt -w 8 --profiling -V --outputMappedCleanReads
```


### Step 3. Database download 
For most non-model organisms, biological understanding of study outcomes is limited to protein-coding genes with functional annotations such as KEGG pathways, Gene Ontology or PANTHER classification system. Therefore, developing Seq2Fun database to focus on functionally annotated genes such as KOs largely meets the preferred needs of most scientists studying non-model organisms.

We provide dozens of pre-built databases that can be downloaded here.
Note: * these KOs in the database are KOs assigned to KEGG pathways and they are only a proportion of whole list of KOs.
** all KOs include KOs not assigned to KEGG pathways.

| Group | Species	| Proteins*	| KOs*	| Filename*	| Proteins**	| KOs**	| Filename** |
| ----- | ------- | --------- | ----- | --------- | ----------- | ----- | ---------- |
| Eukaryotes	| 537	| 1,908,542	| 8,041	| [eukaryotes.tar.gz](http://gofile.me/4esAc/lwO0sfr5c)	| 3,950,549	| 15,302	| [eukaryotes_all_KOs.tar.gz](http://gofile.me/4esAc/3GNJ1lfNg)| 
| Animals	| 250	| 1,126,598	| 6,723	| [animals.tar.gz](http://gofile.me/4esAc/MfQbpPK02)	| 2,446,258	| 12,984	| [animals_all_KOs.tar.gz](http://gofile.me/4esAc/XkLMHAOnn) | 
| Plants	| 105	| 480,379	| 3,012	| [plants.tar.gz](http://gofile.me/4esAc/jp5UT8xwG")	| 926,166	| 6,363	| [plants_all_KOs.tar.gz](http://gofile.me/4esAc/2UTNtmlWs) | 
| Fungi	| 130	| 237,631	| 2,423	| [fungi.tar.gz](http://gofile.me/4esAc/1vo0BgUTU)	| 444,690	| 4,987	| [fungi_all_KOs.tar.gz](http://gofile.me/4esAc/gwqs6EXU0) 
| Protists	| 51	| 64,058	| 2,696	| [protists.tar.gz](http://gofile.me/4esAc/Cbb9LcFLS)	| 133,614	| 6,505	| [protists_all_KOs.tar.gz](http://gofile.me/4esAc/So4gxKC7Q) | 
| Mammals	| 66	| 378,311	| 5,622	| [mammals.tar.gz](http://gofile.me/4esAc/PYOHX6r1y)	| 689,252	| 11,078	| [mammals_all_KOs.tar.gz](http://gofile.me/4esAc/4SGD2Psmp) | 
| Birds	| 24	| 87,530	| 4,177	| [birds.tar.gz](http://gofile.me/4esAc/G0QBrqmvn)	| 208,153	| 9,718	| [birds_all_KOs.tar.gz](http://gofile.me/4esAc/8FtrhiJf0) |
| Reptiles	| 12	| 62,677	| 4,342	| [reptiles.tar.gz](http://gofile.me/4esAc/GjjKjzAsL)	| 153,373	| 10,113	| [reptiles_all_KOs.tar.gz](http://gofile.me/4esAc/K5YUqVVtW) | 
| Amphibians	| 3	| 20,880	| 4,207	| [amphibians.tar.gz](http://gofile.me/4esAc/zUegm99hY)	| 50,137	| 9,715	| [amphibians_all_KOs.tar.gz](http://gofile.me/4esAc/qAmasYGsg) | 
| Fishes	| 39	| 273,691	| 4,308	| [fishes.tar.gz](http://gofile.me/4esAc/VCC2q7LC6)	| 783,801	| 10,510	| [fishes_all_KOs.tar.gz](http://gofile.me/4esAc/Iwjwuxoym) | 
| Arthropods	| 72	| 196,277	| 3,541	| [arthropods.tar.gz](http://gofile.me/4esAc/tu3lbSwk2)	| 455,750	| 8,723	| [arthropods_all_KOs.tar.gz](http://gofile.me/4esAc/EGYQpJobd) | 
| Nematodes	| 6	| 13,379	| 2,324	| [nematodes.tar.gz](http://gofile.me/4esAc/6axLombcd)	| 30,128	| 5,260	| [nematodes_all_KOs.tar.gz](http://gofile.me/4esAc/rS4CkbRpt) | 

## Case Studies
This short tutorial below demonstrates how to run Seq2Fun. We use a RNA-seq dataset from a real non-model organism double-crested cormorant (DCCO), treated with ethinyl estradiol (EE2) as a show case.

### Experimental design
#### 1. Description of experiment

DCCO embryos were exposed via egg injection to EE2, a synthetic estrogen that is the active substance in some forms of birth control. Livers were harvested after 14 days exposure and immediately frozen in liquid nitrogen for total RNA extraction. Total RNA was sent to Genome Quebec (Montreal, Quebec, Canada), to build sequencing library with TruSeq RNA Library Prep Kit (San Diego, California, United States) before submitted to Illumina NovaSeq 6000 (San Diego, California, United States) for 100 bp PE reads sequencing.

#### 2. Experimental samples

Each sample was subsampled with 5 million reads, just for demonstration purpose.
The samples can be download from here

| Group | Chemicals (dose)	| Number of samples	| Number of reads	|
| ----- | ----------------- | ----------------- | --------------- |
| High	| EE2 (33 mg/ml)	  |    4	            |     20,000,000	|
| Medium| EE2 (3.3 mg/ml)	  |    5	            |     25,000,000	|
| Control| DMSO         	  |    5	            |     25,000,000	|


#### 3. Database used

| Group | Number of proteins	| Number of KOs	| Number of species	| Database name |
| ----- | ----------------    | ------------- | ---------------   | ------------- |
| Birds	| 87,530          	  |    4,177	    |     24         	  |[birds.tar.gz](http://gofile.me/4esAc/G0QBrqmvn)|   

### Running Seq2Fun
#### 1.Preparing sample.txt file

This file consists of 3 columns and separated by '\t'. The first column is the prefix name of each sample, and the second is the forward reads file, the last column is the reverse reads file. If you have a single-end (SE) reads, remove the reverse reads (second) column. It looks like this:

A1.CE2-S1-LT	A1.CE2-S1-LT_R1.fastq.gz	A1.CE2-S1-LT_R2.fastq.gz

A2.CE2-M4-LT	A2.CE2-M4-LT_R1.fastq.gz	A2.CE2-M4-LT_R2.fastq.gz

B1.CE2-S2-LT	B1.CE2-S2-LT_R1.fastq.gz	B1.CE2-S2-LT_R2.fastq.gz

B2.CE2-M5-LT	B2.CE2-M5-LT_R1.fastq.gz	B2.CE2-M5-LT_R2.fastq.gz

C1.CE2-S3-LT	C1.CE2-S3-LT_R1.fastq.gz	C1.CE2-S3-LT_R2.fastq.gz

D1.CE2-S4-LT	D1.CE2-S4-LT_R1.fastq.gz	D1.CE2-S4-LT_R2.fastq.gz

D2.CE2-H2-LT	D2.CE2-H2-LT_R1.fastq.gz	D2.CE2-H2-LT_R2.fastq.gz

E1.CE2-S5-LT	E1.CE2-S5-LT_R1.fastq.gz	E1.CE2-S5-LT_R2.fastq.gz

E2.CE2-H3-LT	E2.CE2-H3-LT_R1.fastq.gz	E2.CE2-H3-LT_R2.fastq.gz

F1.CE2-M1-LT	F1.CE2-M1-LT_R1.fastq.gz	F1.CE2-M1-LT_R2.fastq.gz

F2.CE2-H4-LT	F2.CE2-H4-LT_R1.fastq.gz	F2.CE2-H4-LT_R2.fastq.gz

G1.CE2-M2-LT	G1.CE2-M2-LT_R1.fastq.gz	G1.CE2-M2-LT_R2.fastq.gz

G2.CE2-H5-LT	G2.CE2-H5-LT_R1.fastq.gz	G2.CE2-H5-LT_R2.fastq.gz

H1.CE2-M3-LT	H1.CE2-M3-LT_R1.fastq.gz	H1.CE2-M3-LT_R2.fastq.gz






## Bugs or feature requests

To inform us of any bugs or requests, please open a new issue or send an email to #jasmine.chong@mail.mcgill.ca.

## MicrobiomeAnalystR History & Updates

11-16-2020 - Code update w. web + change files from .rds to .qs - users need to install qs R package now
02-24-2020 - Code update w. web + added note about usage
09-05-2019 - Bug fixing w. web
08-07-2019 - Added function to import SILVA annotated biom files (handling Domain in taxonomy)
07-11-2019 - Added volcano + dot plots for RNAseq analysis
07-08-2019 - Testing R code for local use + creating vignettes
07-03-2019 - Updating R code + documentation
06-22-2019 - Prepping R package for stable release
