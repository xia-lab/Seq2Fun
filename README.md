## Description 

**Seq2Fun version 2** is an ultra-fast, all-in-one high resolution functional profiling tool directly from RNA-seq raw reads for organisms without reference genomes. For more detailed descriptions of the concept and algorithm, and instructions, please visit the Seq2Fun website **www.seq2fun.ca**.

If you want to use **version 1** , please download the source code [here](http://gofile.me/4esAc/OY02RFz5y) and database from [here](http://gofile.me/4esAc/UPNLXPxxO)


## Key features of Seq2Fun version 2

* **Ultra-fast**: > 120 times faster (~ 2 million reads / minute) than the conventional RNA-seq workflow.

* **Extremely low memory cost**: consumes as little as 2.27 GB memory and can run on a standard PC with 8 threads and 16 GB memory.

* **Reference-free**: does not require the genome or transcriptome reference of the organism; it is also transcriptome de novo assembly-free.

* **All-in-one**: Seq2Fun directly takes raw RNA-seq reads as input and output gene abundance table without any intermediate file writing and loading, making I/O very efficient.

* **Multifunctional**: Seq2Fun generates several output files, including ortholog gene bundance tables, a html report summarizing these tables and reads quality check, as well as output mapped clean reads for further analysis such as gene assembly.

* **Flexible**: supports RNA-seq analysis on particular genes or groups of organisms using customized database.

* **Easy to use**: requires minimal programing skills.
* **Support target gene assemble**: new features to extract mapped reads and conduct target gene assemble.

## Getting started
### Step 1. Install the software
Seq2Fun (version 2.0.0) is written in C/C++11 and can be installed on Linux or Mac OS X (with Xcode and Xcode Command Line Tools installed). 
We have tested Seq2Fun on Ubuntu (16.04 LTS and above) and macOS Catalina.

```
git clone https://github.com/xia-lab/Seq2Fun.git
cd Seq2Fun/src/
make clean
make
```
### Step 2. Database download 
For most non-model organisms, biological understanding of study outcomes is limited to protein-coding genes with functional annotations such as KEGG pathways, Gene Ontology or PANTHER classification system. Therefore, developing Seq2Fun version 2 database to focus on functionally annotated genes such as orthologs largely meets the preferred needs of most scientists studying non-model organisms.

We provide dozens of pre-built databases that can be downloaded here based on [Orthofinder](https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_orthologs.gz).

| Group | Species | Proteins	| n. of ortholog | n. of core ortholog | Filename	|
| ----- | ------- | --------- | ----- | --------- |
| Vertebrates| 212 | 4,573,967	| 74,321 | [vertebrates.tar.gz](http://gofile.me/4esAc/YkygX7SwI)	|
|Algae | 14 | 155495| 34628| 1831| [algae.tar.gz](http://gofile.me/4esAc/YkygX7SwI)|
|alveolates | 21 | 207674| 48656| 1099| [alveolates.tar.gz](http://gofile.me/4esAc/YkygX7SwI)|
|amoebozoa | 7 |81844| 20217| 1200| [amoebozoa.tar.gz](http://gofile.me/4esAc/YkygX7SwI)|
| amphibians|  3| 75261| 9838| 7776| [amphibians.tar.gz](http://gofile.me/4esAc/YkygX7SwI)|
| animals|  370| 7150735| 236447| 3388| [animals.tar.gz](http://gofile.me/4esAc/YkygX7SwI)|
| apicomplexans|  18|93576| 13823| 1042| [apicomplexans.tar.gz](http://gofile.me/4esAc/YkygX7SwI)|
| arthropods|  119| 1727651| 101867| 4053| [arthropods.tar.gz](http://gofile.me/4esAc/YkygX7SwI)|
| ascomycetes| 100 |904642| 93433| 2590| [ascomycetes.tar.gz](http://gofile.me/4esAc/YkygX7SwI)|
| basidiomycetes| 33 |363997| 52723| 2263| [basidiomycetes.tar.gz](http://gofile.me/4esAc/YkygX7SwI)|
| birds|  31|482205|13868 | 6433| [birds.tar.gz](http://gofile.me/4esAc/YkygX7SwI)|
| cnidarians| 9 |203000| 18682|4840 | [cnidarians.tar.gz](http://gofile.me/4esAc/YkygX7SwI)|
| crustaceans|7|154960|32225 | 4408*| [crustaceans.tar.gz](http://gofile.me/4esAc/YkygX7SwI)|
| dothideomycetes|  10|123200| 26288|4644 | [dothideomycetes.tar.gz](http://gofile.me/4esAc/YkygX7SwI)|
| eudicots|93| 3180221	|72086 | 6935| [eudicots.tar.gz](http://gofile.me/4esAc/YkygX7SwI)|
| euglenozoa|9|86483| 11790| 4341| [euglenozoa.tar.gz](http://gofile.me/4esAc/YkygX7SwI)|
| eurotiomycetes| 20 | 196228	| 23006| 4081| [eurotiomycetes.tar.gz](http://gofile.me/4esAc/YkygX7SwI)|
| fishes|64|1736572| 27327| 7140| [fishes.tar.gz](http://gofile.me/4esAc/YkygX7SwI)|
| flatworms|4|58181| 15156| 3720| [flatworms.tar.gz](http://gofile.me/4esAc/YkygX7SwI)|
|fungi |138|1278312| 141223| 2063| [fungi.tar.gz](http://gofile.me/4esAc/YkygX7SwI)|
| insects| 101 |1376824| 60375| 4508| [insects.tar.gz](http://gofile.me/4esAc/YkygX7SwI)|
| leotiomycetes|  5|67865|19335 | 4833| [leotiomycetes.tar.gz](http://gofile.me/4esAc/YkygX7SwI)|
| mammals|94|1910363| 32471| 7781*| [mammals.tar.gz](http://gofile.me/4esAc/YkygX7SwI)|
| mollusks|9|206905|29992|5485** | [mollusks.tar.gz](http://gofile.me/4esAc/YkygX7SwI)|
| monocots|17|560027|32745| 6060| [monocots.tar.gz](http://gofile.me/4esAc/YkygX7SwI)|
|nematodes | 6 |134093| 32549| 3018| [nematodes.tar.gz](http://gofile.me/4esAc/YkygX7SwI)|
| plants|127|3968027|128122 | 3872| [plants.tar.gz](http://gofile.me/4esAc/YkygX7SwI)|
| protists|52|660237| 126969| 736| [protists.tar.gz](http://gofile.me/4esAc/YkygX7SwI)|
| reptiles| 20 |384584| 11946|6786 | [reptiles.tar.gz](http://gofile.me/4esAc/YkygX7SwI)|
| saccharomycetes| 36 |195913|13337 | 2741| [saccharomycetes.tar.gz](http://gofile.me/4esAc/YkygX7SwI)|
| stramenopiles| 8 |119746| 28902| 680| [stramenopiles.tar.gz](http://gofile.me/4esAc/YkygX7SwI)|
| vertebrates|212|4588985| 59392| 6045| [vertebrates.tar.gz](http://gofile.me/4esAc/YkygX7SwI)|



### Step 3. Run a small test 
There are four sub folders under seq2fun - src, bin, database and testdata. The bin folder contains the binary code we just complied. The testdata contains a small test data from the Case Study.

Download the [birds database](http://gofile.me/4esAc/bG1RlcQIS) to the <code><b> testdata </b></code> folder, and issue the following commands:
```
tar -xzvf birds.tar.gz
```

From within the <code><b> testdata </b></code> folder, issue the following commands:
```
../bin/seq2fun --sampletable sample.txt --tfmi birds/birds.fmi --genemap birds/birds_annotation.txt -w 8 --profiling -V --outputMappedCleanReads --outputReadsAnnoMap
or if you want the trim the first 6 bases
../bin/seq2fun --sampletable sample.txt --tfmi birds/birds.fmi --genemap birds/birds_annotation.txt --trim_front1 6 --trim_front2 6 -w 8 --profiling -V --outputMappedCleanReads --outputReadsAnnoMap
```

## Case Studies
This short tutorial below demonstrates how to run Seq2Fun. We use a RNA-seq dataset from a real non-model organism double-crested cormorant (DCCO), treated with ethinyl estradiol (EE2) as a show case.

### Experimental design
#### 1. Description of experiment

DCCO embryos were exposed via egg injection to EE2, a synthetic estrogen that is the active substance in some forms of birth control. Livers were harvested after 14 days exposure and immediately frozen in liquid nitrogen for total RNA extraction. Total RNA was sent to Genome Quebec (Montreal, Quebec, Canada), to build sequencing library with TruSeq RNA Library Prep Kit (San Diego, California, United States) before submitted to Illumina NovaSeq 6000 (San Diego, California, United States) for 100 bp PE reads sequencing.

#### 2. Experimental samples

Each sample was subsampled with 5 million reads, just for demonstration purpose.
The samples can be download from [here](http://gofile.me/4esAc/UPNLXPxxO)

| Group | Chemicals (dose)	| Number of samples	| Number of reads	|
| ----- | ----------------- | ----------------- | --------------- |
| High	| EE2 (31.9 mg/ml)	  |    5	            |     25,000,000	|
| Medium| EE2 (2.3 mg/ml)	  |    5	            |     25,000,000	|
| Control| DMSO         	  |    5	            |     25,000,000	|


#### 3. Database used

| Group | Number of species | Number of proteins | Number of ortholog | Database name |
| ----- | ----------------    | ------------- | ---------------   | ------------- |
| Birds	| 31 | 482,205 | 24,076 | [birds.tar.gz](http://gofile.me/4esAc/YkygX7SwI) |

### Running Seq2Fun
#### 1.Preparing sample.txt file

This file consists of 4 columns and separated by '\t'. The first column is the prefix name of each sample, and the second is the forward reads file, the thrid column is the reverse reads file, and last one is the sample group info. If you have a single-end (SE) reads, remove the reverse reads (third) column. It looks like this:
```
A1.CE2-S1	A1.CE2-S1_R1.fastq.gz	A1.CE2-S1_R1.fastq.gz	control
A3.CE1-M2	A3.CE1-M2_R1.fastq.gz	A3.CE1-M2_R1.fastq.gz	middle
A4.CE1-H5	A4.CE1-H5_R1.fastq.gz	A4.CE1-H5_R1.fastq.gz	high
B1.CE2-S2	B1.CE2-S2_R1.fastq.gz	B1.CE2-S2_R1.fastq.gz	control
B3.CE1-M3	B3.CE1-M3_R1.fastq.gz	B3.CE1-M3_R1.fastq.gz	middle
C1.CE2-S3	C1.CE2-S3_R1.fastq.gz	C1.CE2-S3_R1.fastq.gz	control
C3.CE1-M4	C3.CE1-M4_R1.fastq.gz	C3.CE1-M4_R1.fastq.gz	middle
D1.CE2-S4	D1.CE2-S4_R1.fastq.gz	D1.CE2-S4_R1.fastq.gz	control
D3.CE1-M5	D3.CE1-M5_R1.fastq.gz	D3.CE1-M5_R1.fastq.gz	middle
E1.CE2-S5	E1.CE2-S5_R1.fastq.gz	E1.CE2-S5_R1.fastq.gz	control
E3.CE1-H1	E3.CE1-H1_R1.fastq.gz	E3.CE1-H1_R1.fastq.gz	high
F3.CE1-H2	F3.CE1-H2_R1.fastq.gz	F3.CE1-H2_R1.fastq.gz	high
G3.CE1-H3	G3.CE1-H3_R1.fastq.gz	G3.CE1-H3_R1.fastq.gz	high
H2.CE1-M1	H2.CE1-M1_R1.fastq.gz	H2.CE1-M1_R1.fastq.gz	middle
H3.CE1-H4	H3.CE1-H4_R1.fastq.gz	H3.CE1-H4_R1.fastq.gz	high
```

#### 2. Running Seq2Fun to quantify RNA-seq reads.

Seq2Fun has two output modes: comparative and profiling mode (default).
The comparative mode produces only the ortholog abundance table, while the profiling mode produces 4 tables 1). ortholog abundance table for all samples, 2). ortholog abundance table for all samples submit to networkanalyst and 3). annotation file sumbit to networkanalyst for downstream analysis, and ortholog abundance table for each sample, 4). mapped clean reads file, and 5). a html report summarizing these tables.
```
S2F_HOME/bin/seq2fun --sampletable sample.txt --tfmi S2F_HOME/database/birds/birds.fmi --genemap S2F_HOME/database/birds/birds_annotation.txt -w 8 --profiling --outputMappedCleanReads
of if you want to trim the first 6 bases
S2F_HOME/bin/seq2fun --sampletable sample.txt --tfmi S2F_HOME/database/birds/birds.fmi --genemap S2F_HOME/database/birds/birds_annotation.txt --trim_front1 6 --trim_front2 6 -w 8 --profiling --outputMappedCleanReads
```

### Main results
#### 1  abundance for all the samples table (S2fid_abundance_table_all_samples.txt)

This table has sample names, meta data and s2fid, reads count (how many reads have assigned to the s2fid/ortholog), KO, GO, gene symbol and gene description separated by '\t'.
```
#NAME	A1.CE2-S1	A3.CE1-M2	A4.CE1-H5	B1.CE2-S2	B3.CE1-M3	C1.CE2-S3	C3.CE1-M4	D1.CE2-S4	D3.CE1-M5	E1.CE2-S5	E3.CE1-H1	F3.CE1-H2	G3.CE1-H3	H2.CE1-M1	H3.CE1-H4	annotation
#CLASS:XX	control	middle	high	control	middle	control	middle	control	middle	control	high	high	high	middle	high	-
s2f_10	596	723	689	326	721	728	831	459	407	394	616	633	499	823	681	K13524|GO:0003824;GO:0003867;GO:0005739;GO:0008483;GO:0009448;GO:0009450;GO:0016740;GO:0030170;GO:0032144;GO:0034386;GO:0042135;GO:0042802;GO:0047298;GO:0048148;|ABAT|4-aminobutyrate aminotransferase
s2f_100	6	17	22	1	22	7	27	5	10	6	13	11	5	17	24	K04137|GO:0004930;GO:0004935;GO:0004937;GO:0005886;GO:0007165;GO:0007186;GO:0016020;GO:0016021;GO:0071875;|ADRA1D|adrenoceptor alpha 1d
s2f_1000	164	107	103	66	148	143	140	104	94	129	92	91	105	123	147	K00344|GO:0008270;GO:0016491|CRYZ|crystallin zeta
s2f_10000	39	52	56	51	77	46	68	45	58	42	72	77	61	79	60	K10105|GO:0006464;GO:0009249|LIPT1|lipoyltransferase 1
s2f_10001	156	323	283	211	324	134	457	214	439	202	368	362	335	309	259	K14565|GO:0001094;GO:0001650;GO:0005634;GO:0005654;GO:0005730;GO:0005732;GO:0005829;GO:0015030;GO:0030515;GO:0031428;GO:0032040;GO:0042254;GO:0048254;GO:0051117;GO:0070761;|NOP58|nop58 ribonucleoprotein
s2f_10002	46	79	83	54	85	46	101	63	122	39	76	90	90	89	48	K25166|GO:0008168;GO:0016740;GO:0032259;|METTL13|methyltransferase 13, eef1a lysine and n-terminal methyltransferase
s2f_10003	53	110	113	126	98	46	110	63	132	45	125	132	105	125	108	K05292|GO:0005737;GO:0016255;GO:0030182;GO:0031410;GO:0042765;GO:0051402;|PIGT|phosphatidylinositol glycan anchor biosynthesis class t
s2f_10004	43	84	77	85	71	57	81	72	117	58	64	93	87	77	82	K03256|GO:0005634;GO:0008033;GO:0030488;GO:0031515;GO:0080009;|TRMT6|trna methyltransferase 6 non-catalytic subunit
s2f_10005	58	92	96	281	83	57	79	63	86	51	124	119	92	93	112	K02144|GO:0000221;GO:0006811;GO:0046961;GO:1902600;|ATP6V1H|atpase h+ transporting v1 subunit h
s2f_10006	123	167	174	149	185	87	167	119	149	72	196	196	143	179	198	K23387|GO:0045048|GET4|guided entry of tail-anchored proteins factor 4
s2f_10007	28	67	45	35	51	35	52	47	49	34	60	53	74	50	50	K00586|GO:0004164;GO:0008168;GO:0016740;GO:0017183;GO:0032259;|DPH5|diphthamide biosynthesis 5
s2f_10008	125	124	118	149	143	100	177	127	180	91	130	159	147	136	122	K20367|GO:0005783;GO:0006888;GO:0006890;GO:0016020;GO:0016021;GO:0016192;GO:0030134;GO:0030173;GO:0030176;GO:0033116;|ERGIC3|ergic and golgi 3
s2f_10009	14	37	47	46	44	21	52	48	53	26	38	48	39	35	43	K14535|GO:0005634;GO:0006352;GO:0046982|TAF9B|tata-box binding protein associated factor 9b
s2f_1001	110	226	244	674	228	108	223	218	461	159	216	277	353	180	205	K01647|GO:0005759;GO:0016740;GO:0046912|CS|citrate synthase
s2f_10010	1	3	3	34	5	1	2	3	30	1	4	8	8	1	2	U|GO:0007212;GO:0016020;GO:0016021;GO:0032051;GO:0048268;|NSG2|neuronal vesicle trafficking associated 2
s2f_10012	1	1	3	3	2	2	5	0	1	2	1	1	4	2	2	K09208|GO:0000978;GO:0000981;GO:0001228;GO:0003677;GO:0006357;GO:0008285;GO:0045647;GO:0045944;GO:1990837;|KLF13|kruppel like factor 13
s2f_10013	95	138	127	141	138	106	129	93	142	95	173	139	157	121	133	U|GO:0005764;GO:0005765;GO:0016020;GO:0016192;GO:0035658;GO:0043231;|CCZ1|ccz1 homolog, vacuolar protein trafficking and biogenesis associated
s2f_10014	5	11	17	36	13	10	18	14	30	12	17	28	27	20	15	K10417|GO:0005737;GO:0005813;GO:0005815;GO:0005856;GO:0005868;GO:0005874;GO:0005929;GO:0005930;GO:0007368;GO:0031514;GO:0035721;GO:0035735;GO:0035869;GO:0036064;GO:0045177;GO:0045504;GO:0060271;GO:1902017;|DYNC2LI1|dynein cytoplasmic 2 light intermediate chain 1
s2f_10015	197	248	217	197	207	173	233	173	193	168	225	238	218	216	176	K15119|GO:0016020;GO:0016021;GO:0055085;|SLC25A39|solute carrier family 25 member 39
s2f_10017	18	37	55	34	48	24	49	36	67	20	50	44	55	57	48	K18342|GO:0004843;GO:0006508;GO:0008233;GO:0008234|OTUD6B|otu K23387|GO:0045048|GET4|guided entry of tail-anchored proteins fact
...
```
#### 2. abundance file for each sample (eg. A1.CE2-S1_s2fid_abundance.txt).
This table has three columns separated by '\t', s2f_id/ortholog, reads_count(how many reads mapped to the s2f_id) and annotation.
```
#s2f_id	Reads_cout	annotation
s2f_10	596	K13524|GO:0003824;GO:0003867;GO:0005739;GO:0008483;GO:0009448;GO:0009450;GO:0016740;GO:0030170;GO:0032144;GO:0034386;GO:0042135;GO:0042802;GO:0047298;GO:0048148;|ABAT|4-aminobutyrate aminotransferase
s2f_100	6	K04137|GO:0004930;GO:0004935;GO:0004937;GO:0005886;GO:0007165;GO:0007186;GO:0016020;GO:0016021;GO:0071875;|ADRA1D|adrenoceptor alpha 1d
s2f_1000	164	K00344|GO:0008270;GO:0016491|CRYZ|crystallin zeta
s2f_10000	39	K10105|GO:0006464;GO:0009249|LIPT1|lipoyltransferase 1
s2f_10001	156	K14565|GO:0001094;GO:0001650;GO:0005634;GO:0005654;GO:0005730;GO:0005732;GO:0005829;GO:0015030;GO:0030515;GO:0031428;GO:0032040;GO:0042254;GO:0048254;GO:0051117;GO:0070761;|NOP58|nop58 ribonucleoprotein
s2f_10002	46	K25166|GO:0008168;GO:0016740;GO:0032259;|METTL13|methyltransferase 13, eef1a lysine and n-terminal methyltransferase
s2f_10003	53	K05292|GO:0005737;GO:0016255;GO:0030182;GO:0031410;GO:0042765;GO:0051402;|PIGT|phosphatidylinositol glycan anchor biosynthesis class t
s2f_10004	43	K03256|GO:0005634;GO:0008033;GO:0030488;GO:0031515;GO:0080009;|TRMT6|trna methyltransferase 6 non-catalytic subunit
s2f_10005	58	K02144|GO:0000221;GO:0006811;GO:0046961;GO:1902600;|ATP6V1H|atpase h+ transporting v1 subunit h
s2f_10006	123	K23387|GO:0045048|GET4|guided entry of tail-anchored proteins factor 4
...                                     ...     ...
```
#### 3. Reads annotation file (eg. A1.CE2-S1_mapped_R1.fastq.gz).
each read was tagged with mapped s2f_id/ortholog
```
@A00266:275:HLFTWDSXX:2:2562:20925:30405 1:N:0:TCGGATTC+CGCAACTA        s2f_3808
CCCGGATAATGAGTGGAACGGAGCAGATGGTGAAGAGGACGGTCATGAGCCCCAGCAGGATGAGGTGGTCGAGCTCCTCCATGCGAGGGGCAGCGGTGGCA
+
FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
@A00266:275:HLFTWDSXX:2:1647:32579:6057 1:N:0:TCGGATTC+CGCAACTA s2f_2261
GTCTGTGACATCCTTTGCCTCAATGTGGCTAGGGGCATCAAGCTTTGTGGTTATCACTCTGGATAGTCCTGGTCCTCTGGTATTGTTTTTCACAATATGAA
+
FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFF:FFFFF
@A00266:275:HLFTWDSXX:2:2608:13711:33098 1:N:0:TCGGATTC+CGCAACTA        s2f_141
CCAGGAGGAAGCCGTTCTTGCAGGGAGGGGAGTCGAGGTTGTTCTCAATCAGCTCCACCACCATCTCATCACTCACCAGCTTCCCCGCATCCATCGTCTCC
+
FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FF,FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
@A00266:275:HLFTWDSXX:2:1327:17996:28714 1:N:0:TCGGATTC+CGCAACTA        s2f_1351
AGCGCAGCACCGAGTAGCTGTTCTCATCGACCTTCTTGCTGTAGGTGACCTCGTACTTCCAGACCCGGCTCTGCTGCCGCGGTGGCACCGTCCAGGAGAGG
+
FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFF
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
      --outputReadsAnnoMap            enable output mapped clean reads-annotation map into .gz files, by default is false, using --outputReadsAnnoMap to enable it

  // Homology search;

   -D, --genemap                    gene/protein species map

       --profiling                  by default it is off. If this option is specified, several output files will be generated, s2f_id/ortholog abundance table.

    
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
12-21-2021 - seq2fun_v2.0.0 released  
08-23-2021 - seq2fun_v1.2.4 released  
06-18-2021 - seq2fun_v1.2.3 released  
06-05-2021 - seq2fun_v1.2.2 released  
03-31-2021 - seq2fun_v1.2.1 released  
03-26-2021 - seq2fun_v1.2.0 released  
03-22-2021 - seq2fun_v1.1.4 released  
01-14-2021 - seq2fun_v1.1.2 released  
12-06-2020 - seq2fun_v1.1.0 released  
08-24-2020 - seq2fun_v1.0.0 released  
