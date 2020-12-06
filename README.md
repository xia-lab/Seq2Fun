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

#### Option B) Clone Github and install locally

The * must be replaced by what is actually downloaded and built. For instance, check your Downloads folder to see what tar.gz file was downloaded. So if you download MicrobiomeAnalystR_1.0.1.tar.gz, replace the * with the downloaded version number.  

```R
git clone https://github.com/xia-lab/MicrobiomeAnalystR.git
R CMD build MicrobiomeAnalystR
R CMD INSTALL MicrobiomeAnalystR_*.tar.gz

```

#### Option C) Manual download of MicrobiomeAnalystR.tar.gz and install locally - not yet available, stable release to come soon!!

Manually download the .tar.gz file from [here](https://www.dropbox.com/s/wk43rs9hswzypgt/MicrobiomeAnalystR_0.0.0.9000.tar.gz?dl=0). The * must be replaced by what is actually downloaded and built.  

```R
cd ~/Downloads
R CMD INSTALL MicrobiomeAnalystR_*.tar.gz

```

## Case Studies

### MicrobiomeAnalyst Workflow

To showcase how to utilize MicrobiomeAnalystR , we provide a detailed tutorial to perform a comprehensive end-to-end workflow from raw sequence data preprocessing to knowledge-based analysis. The dataset showcased in the tutorial consists of a subset of pediatric IBD stool samples obtained from the Integrative Human Microbiome Project Consortium (https://ibdmdb.org/). The tutorial is available inside the R package as a vignette.

## Tutorials

For detailed tutorials on how to use MicrobiomeAnalystR, please refer to the R package vignettes. These vignettes include a comprehensive tutorial introducing MicrobiomeAnalystR, four detailed step-by-step tutorials with example data for each of the main MetaboAnalytR  modules, and a case-study showcasing the end-to-end functionality of MicrobiomeAnalystR. Note, the functions below work only if the R package vignettes were built. 

Within R:
```R
vignette(package="MicrobiomeAnalystR")
```

Within a web-browser:
```R
browseVignettes("MicrobiomeAnalystR")
```

## Citation

MicrobiomeAnalystR has been developed by the [XiaLab](http://xialabresearch.com/) at McGill University. The original manuscript (web-based version) can be found [here](https://www.ncbi.nlm.nih.gov/pubmed/28449106). 

We encourage users to further develop the package to suit their needs. If you use the R package, please cite us: 

Dhariwal A, Chong J, Habib S, King IL, Agellon LB, Xia J. MicrobiomeAnalyst: a web-based tool for comprehensive statistical, visual and meta-analysis of microbiome data. Nucleic acids research. 2017 Jul 3;45(W1):W180-8.

*From within R:*

```R
citation("MicrobiomeAnalystR")
```

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

## MicrobiomeAnalystR TO DO


