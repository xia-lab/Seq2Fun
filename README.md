## Description 

**Seq2Fun version 2** is an ultra-fast, all-in-one functional profiling tool directly from RNA-seq raw reads for organisms without reference genomes. For more detailed descriptions of the concept and algorithm, and instructions, please visit the Seq2Fun website **www.seq2fun.ca** .


## Key features of Seq2Fun version 2

* **Ultra-fast**: Seq2Fun is > 120 times faster (~ 2 million reads / minute) than the conventional RNA-seq workflow.

* **Extremely low memory cost**: Seq2Fun consumes as little as 2.27 GB memory and can run on a standard PC with 8 threads and 16 GB memory.

* **Reference-free**: Seq2Fun does not require the genome or transcriptome reference of the organism; it is also transcriptome de novo assembly-free.

* **Highly accurate**: Seq2Fun generates KO abundance with R2 value as high as 0.93 comparing with the ground truth.

* **All-in-one**: Seq2Fun directly takes raw RNA-seq reads as input and output gene abundance table without any intermediate file writing and loading, making I/O very efficient.

* **Multifunctional**: Seq2Fun generates 6 levels of output files, including KO and GO abundance tables, hit pathway table, hit species table, reads KO table, a html report summarizing these tables and reads quality check, as well as output mapped clean reads for further analysis such as gene assembly.

* **Flexible**: Seq2Fun supports RNA-seq analysis on particular genes or groups of organisms using customized database.

* **Easy to use**: Seq2Fun requires minimal programing skills.

## Getting started
### Step 1. Install the package
Seq2Fun (version 2.0.0) is written in C/C++11 and can be installed on Linux or Mac OS X (with Xcode and Xcode Command Line Tools installed). 
We have tested Seq2Fun on Ubuntu (16.04 LTS and above) and macOS Catalina.

```
git clone https://github.com/xia-lab/Seq2Fun.git
cd Seq2Fun/src/
make clean
make
```
### Step 2. Run a small test 
There are four sub folders under seq2fun - src, bin, database and testdata. The bin folder contains the binary code we just complied. The testdata contains a small test data from the Case Study.
From within the <code><b> testdata </b></code> folder, issue the following commands:
```
../bin/seq2fun --sampletable sample.txt --tfmi birds_cdhit99_proteins.fmi --genemap birds_protein_ko_species_cdhit99.txt -w 8 --profiling -V --outputMappedCleanReads --outputReadsKOMap
or if you want the trim the first 6 bases
../bin/seq2fun --sampletable sample.txt --tfmi birds_cdhit99_proteins.fmi --genemap birds_protein_ko_species_cdhit99.txt --trim_front1 6 --trim_front2 6 -w 8 --profiling -V --outputMappedCleanReads --outputReadsKOMap
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
| High	| EE2 (31.9 mg/ml)	  |    5	            |     25,000,000	|
| Medium| EE2 (2.3 mg/ml)	  |    5	            |     25,000,000	|
| Control| DMSO         	  |    5	            |     25,000,000	|


#### 3. Database used

| Group | Number of proteins	| Number of KOs	| Number of species	| Database name |
| ----- | ----------------    | ------------- | ---------------   | ------------- |
| Birds	| 87,530          	  |    4,177	    |     24         	  |[birds.tar.gz](http://gofile.me/4esAc/G0QBrqmvn)|   

### Running Seq2Fun
#### 1.Preparing sample.txt file

This file consists of 4 columns and separated by '\t'. The first column is the prefix name of each sample, and the second is the forward reads file, the thrid column is the reverse reads file, and last one is the sample group info. If you have a single-end (SE) reads, remove the reverse reads (third) column. It looks like this:
```
A1.CE2-S1-LT	A1.CE2-S1-LT_R1.fastq.gz	A1.CE2-S1-LT_R2.fastq.gz	control
A2.CE2-M4-LT	A2.CE2-M4-LT_R1.fastq.gz	A2.CE2-M4-LT_R2.fastq.gz	medium
B1.CE2-S2-LT	B1.CE2-S2-LT_R1.fastq.gz	B1.CE2-S2-LT_R2.fastq.gz	control
B2.CE2-M5-LT	B2.CE2-M5-LT_R1.fastq.gz	B2.CE2-M5-LT_R2.fastq.gz	medium
C1.CE2-S3-LT	C1.CE2-S3-LT_R1.fastq.gz	C1.CE2-S3-LT_R2.fastq.gz	control
D1.CE2-S4-LT	D1.CE2-S4-LT_R1.fastq.gz	D1.CE2-S4-LT_R2.fastq.gz	control
D2.CE2-H2-LT	D2.CE2-H2-LT_R1.fastq.gz	D2.CE2-H2-LT_R2.fastq.gz	high
E1.CE2-S5-LT	E1.CE2-S5-LT_R1.fastq.gz	E1.CE2-S5-LT_R2.fastq.gz	control
E2.CE2-H3-LT	E2.CE2-H3-LT_R1.fastq.gz	E2.CE2-H3-LT_R2.fastq.gz	high
F1.CE2-M1-LT	F1.CE2-M1-LT_R1.fastq.gz	F1.CE2-M1-LT_R2.fastq.gz	medium
F2.CE2-H4-LT	F2.CE2-H4-LT_R1.fastq.gz	F2.CE2-H4-LT_R2.fastq.gz	high
G1.CE2-M2-LT	G1.CE2-M2-LT_R1.fastq.gz	G1.CE2-M2-LT_R2.fastq.gz	medium
G2.CE2-H5-LT	G2.CE2-H5-LT_R1.fastq.gz	G2.CE2-H5-LT_R2.fastq.gz	high
H1.CE2-M3-LT	H1.CE2-M3-LT_R1.fastq.gz	H1.CE2-M3-LT_R2.fastq.gz	medium
```

#### 2. Running Seq2Fun to quantify RNA-seq reads.

Seq2Fun has two output modes: comparative and profiling mode (default).
The comparative mode produces only the KO abundance table, while the profiling mode produces 4 tables 1). KO abundance table for all samples, and KO abundance table for each sample, 2). hit pathway, 3). hit species, 4). reads KO table, and 5). a html report summarizing these tables summarizing these tables.
```
S2F_HOME/bin/seq2fun --sampletable sample.txt --tfmi S2F_HOME/database/birds/birds_cdhit99_proteins.fmi --genemap S2F_HOME/database/birds/birds_protein_ko_species_cdhit99.txt -w 8 --profiling
of if you want to trim the first 6 bases
S2F_HOME/bin/seq2fun --sampletable sample.txt --tfmi S2F_HOME/database/birds/birds_cdhit99_proteins.fmi --genemap S2F_HOME/database/birds/birds_protein_ko_species_cdhit99.txt --trim_front1 6 --trim_front2 6 -w 8 --profiling
```

### Results
#### 1 KO abundance for all the samples table (KO_abundance.txt)

This table has KO id, sample names and KO name separated by '\t'. (how many reads have assigned to the homology KO), the full name of the assigned KO. It looks like this:
```
#Name   A1.CE2-S1-LT    A2.CE2-M4-LT    B1.CE2-S2-LT    B2.CE2-M5-LT    C1.CE2-S3-LT    D1.CE2-S4-LT    D2.CE2-H2-LT    E1.CE2-S5-LT    E2.CE2-H3-LT    F1.CE2-M1-LT    F2.CE2-H4-LT    G1.CE2-M2-LT    G2.CE2-H5-LT    H1.CE2-M3-LT    KO_name
#Class	control		medium		control		medium		control		control		high		control		high		medium		high		medium		high		medium		class_info	
K00002	118             96              386             131             147             141             106             129             120             98              148             117             136             121             AKR1A1, adh; alcohol dehydrogenase (NADP+) [EC:1.1.1.2]
K00006	629             604             235             648             664             506             628             499             670             455             838             615             579             521             GPD1; glycerol-3-phosphate dehydrogenase (NAD+) [EC:1.1.1.8]
K00008	971             755             715             770             1122            770             1058            1010            1023            829             1055            1351            954             1139            SORD, gutB; L-iditol 2-dehydrogenase [EC:1.1.1.14]
K00010	17              18              31              17              17              14              15              17              17              17              15              9               20              25              iolG; myo-inositol 2-dehydrogenase   D-chiro-inositol 1-dehydrogenase [EC:1.1.1.18 1.1.1.369]
K00011	276             292             1581            303             315             336             290             305             353             263             316             279             290             296             AKR1B; aldehyde reductase [EC:1.1.1.21]
K00012	130             94              609             112             111             134             127             344             103             193             163             382             116             299             UGDH, ugd; UDPglucose 6-dehydrogenase [EC:1.1.1.22]
...
```
#### 2. Hit pathway table (A1.CE2-S1_pathway_hits.txt).

This table has five columns separated by '\t', pathway_id, pathway_name, KO_id (which hit KOs are mapped to this pathway), KO_count and KO_name. It looks like this:
```
pathway_id	pathway_name                    KO_id	KO_count	KO_name
map00010	Glycolysis_/_Gluconeogenesis	K00002	118             AKR1A1, adh; alcohol dehydrogenase (NADP+) [EC:1.1.1.2]
map00010	Glycolysis_/_Gluconeogenesis	K00016	12019           LDH, ldh; L-lactate dehydrogenase [EC:1.1.1.27]
map00010	Glycolysis_/_Gluconeogenesis	K00121	573             frmA, ADH5, adhC; S-(hydroxymethyl)glutathione dehydrogenase   alcohol dehydrogenase [EC:1.1.1.284 1.1.1.1]
map00010	Glycolysis_/_Gluconeogenesis	K00128	3721            ALDH; aldehyde dehydrogenase (NAD+) [EC:1.2.1.3]
map00010	Glycolysis_/_Gluconeogenesis	K00129	13              E1.2.1.5; aldehyde dehydrogenase (NAD(P)+) [EC:1.2.1.5]
map00010	Glycolysis_/_Gluconeogenesis	K00134	29734           GAPDH, gapA; glyceraldehyde 3-phosphate dehydrogenase [EC:1.2.1.12]
map00010	Glycolysis_/_Gluconeogenesis	K00149	4504            ALDH9A1; aldehyde dehydrogenase family 9 member A1 [EC:1.2.1.47 1.2.1.3]
...             ...                             ...     ...             ...
```
#### 3. Hit species table (A1.CE2-S1_species_hits.txt).

This table has two columns separated by '\t', species name and number of KOs (sorted by descending order) assigned to this species.
```
species                     number_of_KOs
Nipponia_nippon             1470
Egretta_garzetta            1200
Pygoscelis_adeliae          1187
Columba_livia               1057
Athene_cunicularia          1048
Apteryx_mantelli_mantelli   840
Empidonax_traillii          820
Falco_cherrug               779
Gallus_gallus               658
Anas_platyrhynchos          564
Anser_cygnoides_domesticus  542
Falco_peregrinus            510
...                     ...
```
#### 4. Reads KO table (A1.CE2-S1_reads_ko.txt).
This table has three columns separated by '\t', reads_id, KO_id (which homology KO assigned) and KO_name.
```
reads_id                                KO_id	KO_name
A00266:275:HLFTWDSXX:2:1101:10013:26537	K14736	TF; transferrin
A00266:275:HLFTWDSXX:2:1101:10122:16235	K14736	TF; transferrin
A00266:275:HLFTWDSXX:2:1101:10122:17174	K03883	ND5; NADH-ubiquinone oxidoreductase chain 5 [EC:7.1.1.2]
A00266:275:HLFTWDSXX:2:1101:10122:2174	K00602	purH; phosphoribosylaminoimidazolecarboxamide formyltransferase / IMP cyclohydrolase [EC:2.1.2.3 3.5.4.10]
A00266:275:HLFTWDSXX:2:1101:10140:28510	K00799	GST, gst; glutathione S-transferase [EC:2.5.1.18]
A00266:275:HLFTWDSXX:2:1101:10149:34225	K06238	COL6A; collagen, type VI, alpha
A00266:275:HLFTWDSXX:2:1101:10185:19382	K03883	ND5; NADH-ubiquinone oxidoreductase chain 5 [EC:7.1.1.2]
A00266:275:HLFTWDSXX:2:1101:10212:16423	K11188	PRDX6; peroxiredoxin 6, 1-Cys peroxiredoxin [EC:1.11.1.7 1.11.1.15 3.1.1.-]
A00266:275:HLFTWDSXX:2:1101:10212:36777	K08737	MSH6; DNA mismatch repair protein MSH6
A00266:275:HLFTWDSXX:2:1101:10294:9236	K02938	RP-L8e, RPL8; large subunit ribosomal protein L8e
A00266:275:HLFTWDSXX:2:1101:10321:31735	K14736	TF; transferrin
A00266:275:HLFTWDSXX:2:1101:10339:16297	K06171	NCSTN; nicastrin
A00266:275:HLFTWDSXX:2:1101:10375:24972	K04660	DCN; decorin
A00266:275:HLFTWDSXX:2:1101:10402:14121	K02132	ATPeF1A, ATP5A1, ATP1; F-type H+-transporting ATPase subunit alpha
A00266:275:HLFTWDSXX:2:1101:10420:21198	K03231	EEF1A; elongation factor 1-alpha
A00266:275:HLFTWDSXX:2:1101:10438:27179	K00149	ALDH9A1; aldehyde dehydrogenase family 9 member A1 [EC:1.2.1.47 1.2.1.3]
A00266:275:HLFTWDSXX:2:1101:1045:22498	K00255	ACADL; long-chain-acyl-CoA dehydrogenase [EC:1.3.8.8]
A00266:275:HLFTWDSXX:2:1101:1045:33144	K00134	GAPDH, gapA; glyceraldehyde 3-phosphate dehydrogenase [EC:1.2.1.12]
A00266:275:HLFTWDSXX:2:1101:10465:27508	K05692	ACTB_G1; actin beta/gamma 1
...                                     ...     ...
```

## Seq2Fun full usage options
```
Use seq2fun or seq2fun --help to show the full usage options

  options:

  // input/output

  -s, --sampletable,                (recommended) sample table must consist of 3 columns (sample prefix name (sample01), forward reads name (sample01_R1.fq.gz), group info (control) for single-reads or 4 columns (sample prefix name (sample01), forward reads (sample01_R1.fq.gz), reverse reads (sample01_R2.fq.gz), group info (control) for paired-end reads. The columns must be separated by tab 

  -i, --in1                         read1 input file name

  -I, --in2                         read2 input file name

  -X, --prefix                      (not recommended) prefix name for output files, eg: sample01

      --outputMappedCleanReads,          enable output mapped clean reads into fastq.gz files, by default is false, using --outputMappedCleanReads to enable it
      --outputReadsKOMap            enable output mapped clean reads-KO map into .gz files, by default is false, using --outputReadsKOMap to enable it

  // Homology search;

   -D, --genemap                    gene/protein KO species map

       --profiling                  by default it is off. If this option is specified, 4 levels of output files will be generated, ko abundance table, hit pathway table, hit species table and ko reads mapping tableby default it is off. If this option is specified, 4 levels of output files will be generated, ko abundance table, hit pathway table, hit species table and ko reads mapping table

    
  // translated search

   -d, --tfmi                       fmi index of Protein database

   -K, --mode                       searching mode either tGREEDY or tMEM (maximum exactly match). By default greedy 

   -E, --mismatch                   number of mismatched amino acid in sequence comparison with protein database with default value 2

   -j, --minscore                   minimum matching score of amino acid sequence in comparison with protein database with default value 100

   -J, --minlength                  minimum matching length of amino acid sequence in comparison with protein database with default value 25, for GREEDY and MEM model

   -m, --maxtranslength             maximum cutoff of translated peptides, it must be no less than minlength, with default 60
       
       --allFragments               enable this function will force Seq2Fun to use all the translated AA fragments with length > minlength. This will slightly help to classify reads contain the true stop codon and start codon; This could have limited impact on the accuracy for comparative study and enable this function will slow down the Seq2Fun. by default is false, using --allFragments to enable it
       
       --codontable                 select the codon table (same as blastx in NCBI), we provide 21 codon tables from 'https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG31'. By default is the codontable1 (Standard Code), the complete codon table can be seen below.
       
       --dbDir			     specify dir for internal databases such as ko_fullname.txt".
	
    
  //selected pathways

   -Z, --pathway                    list of selected pathways for target pathways analysis

   -z, --genefa                     the gene/protein sequences fasta file for retrieving proteins in selected pathways to construct database

    
  // threading

  -w, --thread                      worker thread number, default is 2

    
  -V, --verbose                     enable verbose

      --debug                       enable debug


      --phred64                     indicate the input is using phred64 scoring (it'll be converted to phred33, so the output will still be phred33)

      --reads_to_process            specify how many reads/pairs to be processed. Default 0 means process all reads


  // adapter

  -A, --disable_adapter_trimming    adapter trimming is enabled by default. If this option is specified, adapter trimming is disabled

  -a, --adapter_sequence            the adapter for read1. For SE data, if not specified, the adapter will be auto-detected. For PE data, this is used if R1/R2 are found not overlapped

      --adapter_sequence_r2         the adapter for read2 (PE data only). This is used if R1/R2 are found not overlapped. If not specified, it will be the same as adapter_sequence

      --adapter_fasta               specify a FASTA file to trim both read1 and read2 (if PE) by all the sequences in this FASTA file

      --detect_adapter_for_pe       by default, the auto-detection for adapter is for SE data input only, turn on this option to enable it for PE data

    
  //polyA tail

      --no_trim_polyA               by default, ployA tail will be trimmed. If this option is specified, polyA trimming is disabled


  // trimming

  -f, --trim_front1                 trimming how many bases in front for read1, default is 0

  -t, --trim_tail1                  trimming how many bases in tail for read1, default is 0

  -b, --max_len1                    if read1 is longer than max_len1, then trim read1 at its tail to make it as long as max_len1. Default 0 means no limitation

  -F, --trim_front2                 trimming how many bases in front for read2. If it's not specified, it will follow read1's settings

  -T, --trim_tail2                  trimming how many bases in tail for read2. If it's not specified, it will follow read1's settings

  -B, --max_len2                    if read2 is longer than max_len2, then trim read2 at its tail to make it as long as max_len2. Default 0 means no limitation. If it's not specified, it will follow read1's settings


  // polyG tail trimming

  -g, --trim_poly_g                 force polyG tail trimming, by default trimming is automatically enabled for Illumina NextSeq/NovaSeq data

      --poly_g_min_len              the minimum length to detect polyG in the read tail. 10 by default

  -G, --disable_trim_poly_g         disable polyG tail trimming, by default trimming is automatically enabled for Illumina NextSeq/NovaSeq data


  // polyX tail trimming

  -x, --trim_poly_x                 enable polyX trimming in 3' ends

      --poly_x_min_len              the minimum length to detect polyX in the read tail. 10 by default


  // cutting by quality

      --cut_front                   move a sliding window from front (5') to tail, drop the bases in the window if its mean quality < threshold, stop otherwise

      --cut_tail                    move a sliding window from tail (3') to front, drop the bases in the window if its mean quality < threshold, stop otherwise

  -r, --cut_right                   move a sliding window from front to tail, if meet one window with mean quality < threshold, drop the bases in the window and the right part, and then stop

  -W, --cut_window_size             the window size option shared by cut_front, cut_tail or cut_sliding. Range: 1~1000, default: 4

  -M, --cut_mean_quality            the mean quality requirement option shared by cut_front, cut_tail or cut_sliding. Range: 1~36 default: 20 (Q20)

      --cut_front_window_size       the window size option of cut_front, default to cut_window_size if not specified

      --cut_front_mean_quality      the mean quality requirement option for cut_front, default to cut_mean_quality if not specified

      --cut_tail_window_size        the window size option of cut_tail, default to cut_window_size if not specified

      --cut_tail_mean_quality       the mean quality requirement option for cut_tail, default to cut_mean_quality if not specified

      --cut_right_window_size       the window size option of cut_right, default to cut_window_size if not specified

      --cut_right_mean_quality      the mean quality requirement option for cut_right, default to cut_mean_quality if not specified


  // quality filtering

  -Q, --disable_quality_filtering   quality filtering is enabled by default. If this option is specified, quality filtering is disabled

  -q, --qualified_quality_phred     the quality value that a base is qualified. Default 15 means phred quality <=Q15 is qualified

  -u, --unqualified_percent_limit   how many percents of bases are allowed to be unqualified (0~100). Default 40 means 40%

  -n, --n_base_limit                if one read's number of N base is >n_base_limit, then this read/pair is discarded. Default is 5

  -e, --average_qual                if one read's average quality score <avg_qual, then this read/pair is discarded. Default 0 means no requirement


  // length filtering

  -L, --disable_length_filtering    length filtering is enabled by default. If this option is specified, length filtering is disabled

  -l, --length_required             reads shorter than length_required will be discarded, default is 60

      --length_limit                reads longer than length_limit will be discarded, default 0 means no limitation


  // low complexity filtering

      --no_low_complexity_filter    disable low complexity filter. The complexity is defined as the percentage of base that is different from its next base (base[i] != base[i+1])

  -Y, --complexity_threshold        the threshold for low complexity filter (0~100). Default is 30, which means 30% complexity is required


  // filter by indexes
      --filter_by_index1            specify a file contains a list of barcodes of index1 to be filtered out, one barcode per line

      --filter_by_index2            specify a file contains a list of barcodes of index2 to be filtered out, one barcode per line

      --filter_by_index_threshold   the allowed difference of index barcode for index filtering, default 0 means completely identical


  // base correction in overlapped regions of paired end data

  -c, --disable_correction          disenable base correction in overlapped regions (only for PE data), default is enabled");

  -v, --overlap_len_require         the minimum length to detect overlapped region of PE reads. This will affect overlap analysis based PE merge, adapter trimming and correction. 30 by default

      --overlap_diff_limit          the maximum number of mismatched bases to detect overlapped region of PE reads. This will affect overlap analysis based PE merge, adapter trimming and correction. 5 by default

      --overlap_diff_percent_limit  the maximum percentage of mismatched bases to detect overlapped region of PE reads. This will affect overlap analysis based PE merge, adapter trimming and correction. Default 20 means 20%


  // umi

  -u, --umi                         enable unique molecular identifier (UMI) preprocessing

      --umi_loc                     specify the location of UMI, can be (index1/index2/read1/read2/per_index/per_read, default is none

      --umi_len                     if the UMI is in read1/read2, its length should be provided

      --umi_prefix                  if specified, an underline will be used to connect prefix and UMI (i.e. prefix=UMI, UMI=AATTCG, final=UMI_AATTCG). No prefix by default

      --umi_skip                    if the UMI is in read1/read2, Seq2Fun can skip several bases following UMI, default is 0


  // overrepresented sequence analysis

  -p, --overrepresentation_analysis enable overrepresented sequence analysis

  -P, --overrepresentation_sampling one in (--overrepresentation_sampling) reads will be computed for overrepresentation analysis (1~10000), smaller is slower, default is 20



  // deprecated options

     --cut_by_quality5              DEPRECATED, use --cut_front instead

     --cut_by_quality3              DEPRECATED, use --cut_tail instead

     --cut_by_quality_aggressive    DEPRECATED, use --cut_right instead

     --discard_unmerged             DEPRECATED, no effect now, see the introduction for merging
```

### codon table
We have followed the codon tables from [NCBI](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi);
| Seq2Fun | NCBI |
| ------- | ---- |
|codontable1 |The Standard Code (transl_table=1)|
|codontable2 |The Vertebrate Mitochondrial Code (transl_table=2)|
|codontable3 |The Yeast Mitochondrial Code (transl_table=3)|
|codontable4 |The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code (transl_table=4)|
|codontable5 |The Invertebrate Mitochondrial Code (transl_table=5)|
|codontable6 |The Ciliate, Dasycladacean and Hexamita Nuclear Code (transl_table=6)|
|codontable9 |The Echinoderm and Flatworm Mitochondrial Code (transl_table=9)|
|codontable10 |The Euplotid Nuclear Code (transl_table=10)|
|codontable12 |The Alternative Yeast Nuclear Code (transl_table=12)|
|codontable13 |The Ascidian Mitochondrial Code (transl_table=13)|
|codontable14 |The Alternative Flatworm Mitochondrial Code (transl_table=14)|
|codontable16 |Chlorophycean Mitochondrial Code (transl_table=16)|
|codontable21 |Trematode Mitochondrial Code (transl_table=21)|
|codontable22 |Scenedesmus obliquus Mitochondrial Code (transl_table=22)|
|codontable24 |Rhabdopleuridae Mitochondrial Code (transl_table=24)|
|codontable26 |Pachysolen tannophilus Nuclear Code (transl_table=26)|
|codontable27 |Karyorelict Nuclear Code (transl_table=27)|
|codontable29 |Mesodinium Nuclear Code (transl_table=29)|
|codontable30 |Peritrich Nuclear Code (transl_table=30)|
|codontable31 |Blastocrithidia Nuclear Code (transl_table=31)|
|codontable33 |Cephalodiscidae Mitochondrial UAA-Tyr Code (transl_table=33)|


## Bugs or feature requests

To inform us of any bugs or requests, please open a new issue or send an email to rocpengliu@gmail.com or jeff.xia@mcgill.ca

## Seq2Fun History & Updates
10-01-2021 - seq2fun_v2.0.0 released  
08-23-2021 - seq2fun_v1.2.4 released  
06-18-2021 - seq2fun_v1.2.3 released  
06-05-2021 - seq2fun_v1.2.2 released  
03-31-2021 - seq2fun_v1.2.1 released  
03-26-2021 - seq2fun_v1.2.0 released  
03-22-2021 - seq2fun_v1.1.4 released  
01-14-2021 - seq2fun_v1.1.2 released  
12-06-2020 - seq2fun_v1.1.0 released  
08-24-2020 - seq2fun_v1.0.0 released  
