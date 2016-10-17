---
title: smRNA analysis notes
layout: post
comments: true
---

I recently analysed some smRNA data for a paper I'm working on. These are my analysis notes.  

I used previously published data for Brachypodium, from this paper:  

Garvin DF, Schmutz J, Rokhsar D, Bevan MW, Barry K, Lucas S, et al. Genome sequencing and analysis of the model grass Brachypodium distachyon. Nature. 2010;463: 763–768. doi:10.1038/nature08747

First step is to download the data:  

```bash
$ wget ftp://ftp-trace.ncbi.nlm.nih.gov//sra/sra-instant/reads/ByStudy/sra/SRP/SRP001/SRP001895/SRR035616/SRR035616.sra
$ fastq-dump SRR035616.sra
$ pigz SRR035616.fastq
```

<!--break-->

It is essential that the reads are trimmed correctly to be able to get useful data out of the reads. To do that, we first need to work out what the adapter sequence is. The easiest way to do that is to run fastqc and see what sequences are overrepresented.

```bash
$ fastqc SRR035616.fastq.gz
```
Then we can do a multiple sequence alignment of overrepresented sequences to find the adapter sequence:  

```bash
CLUSTAL multiple sequence alignment by MUSCLE (3.8)


12              GGTAGTTCGACCGCGGA------------ACTGTAGGCACCATCAATT------------
15              GGTAGTTCGACCGCGGAA-----------TCTGTAGGCACCATCAAT-------------
14              -------TCGCTTGGTGCAGATCGGGA--CCTGTAGGCACCCTCA---------------
18              ---TTCATGGACGTTGATAAGATCCTT--CCTGTATGCACC-------------------
1               -------TCGCTTGGTGCAGATCGGGA--CCTGTAGGCACCATCA---------------
13              -------TGATTGAGCCGTGCCAATAT--CCTGTAGGCACCATCA---------------
3               ------------GGGGGTGTAGCTCAT--ACTGTAGGCACCATCAATTCG----------
10              ------------GGGGATGTAGCTCAG--ACTGTAGGCACCATCAATTCG----------
2               ------------GGGGATGTAGCTCAA--ACTGTAGGCACCATCAATTCG----------
11              --------GACCGCATAGCGCAGTGG---ACTGTAGGCACCATCAAT-------------
9               ---------ATTGAGTGCAGCGTTGATGAACTGTAGGCACCATCA---------------
19              -----ACTGGTTGG--ATCATGCTTCT--ACTGTAGGCACCATCA---------------
8               ----TTTGGATTGAAGGGAGCTCT-----GCTGTAGGCACCATCA---------------
7               ---------TCCACAGGCTTTCTTGAACTGCTGTAGGCACCATCA---------------
17              -----TTGGACTGAAGGGTGCTCCC----TCTGTAGGCACCATCA---------------
4               -----TGAAGCTGCCAGCATGATCTG---ACTGTAGGCACCATC----------------
16              ----ACCTGCTCTGATACCATGTTGTG--ACTGTAGGCACCA------------------
6               ---TTCATGGACGTTGATAAGATCCTT--CCTGTAGGCACC-------------------
5               ------------------------------CTGTAGGCACCATCAATTCGTATGCCGTCT
                                              ***** *****        
```

From this it becomes obvious that the adapter is `CTGTAGGCACCATCA...`. Next step is to trim the adapters, and we will also discard any reads that were not trimmed, as these reads are longer than the longest smRNA (so every read should need trimming). We will also retain a minumum read length of 15 nt.

## Adapter trimming

```
$ cutadapt --discard-untrimmed -m 15 -a CTGTAGGCACCATCAATTCG -o adapertrim.fq SRR035616.fastq.gz

This is cutadapt 1.10 with Python 2.7.10
Command line parameters: --discard-untrimmed -m 15 -a CTGTAGGCACCATCAATTCG -o adapertrim.fq SRR035616.fastq.gz
Trimming 1 adapter with at most 10.0% errors in single-end mode ...
Finished in 68.63 s (15 us/read; 3.94 M reads/minute).

=== Summary ===

Total reads processed:               4,501,847
Reads with adapters:                 3,325,787 (73.9%)
Reads that were too short:              57,000 (1.3%)
Reads written (passing filters):     3,268,787 (72.6%)

Total basepairs processed:   162,066,492 bp
Total written (filtered):     70,887,548 bp (43.7%)

=== Adapter 1 ===

Sequence: CTGTAGGCACCATCAATTCG; Type: regular 3'; Length: 20; Trimmed: 3325787 times.

No. of allowed errors:
0-9 bp: 0; 10-19 bp: 1; 20 bp: 2

Bases preceding removed adapters:
  A: 35.9%
  C: 26.7%
  G: 15.9%
  T: 20.4%
  none/other: 1.0%

Overview of removed sequences
length	count	expect	max.err	error counts
3	1307	70341.4	0	1307
4	503	17585.3	0	503
5	832	4396.3	0	832
6	2505	1099.1	0	2505
7	7373	274.8	0	7373
8	10601	68.7	0	10601
9	23318	17.2	0	20283 3035
10	50136	4.3	1	29746 20390
11	128352	1.1	1	72611 55741
12	1185903	0.3	1	680195 505708
13	272422	0.1	1	168709 103713
14	154255	0.0	1	97083 57172
15	526293	0.0	1	345638 180655
16	140755	0.0	1	84107 56648
17	163358	0.0	1	102567 60791
18	149275	0.0	1	94982 54248 45
19	135400	0.0	1	83441 51271 688
20	251034	0.0	2	120106 89305 41623
21	65165	0.0	2	35913 20238 9014
22	11648	0.0	2	6752 3502 1394
23	3309	0.0	2	2054 876 379
24	2179	0.0	2	1456 518 205
25	1963	0.0	2	1398 416 149
26	1404	0.0	2	983 309 112
27	979	0.0	2	714 193 72
28	768	0.0	2	572 137 59
29	628	0.0	2	496 102 30
30	365	0.0	2	301 52 12
31	192	0.0	2	167 19 6
32	133	0.0	2	97 27 9
33	61	0.0	2	49 10 2
34	39	0.0	2	30 3 6
35	39	0.0	2	29 6 4
36	33293	0.0	2	27424 5411 458
```

Now rerun fastqc

```bash
$ fastqc adapertrim.fq
```

We can see peaks at 21 and 24 nt, which are the most abundant smRNA size classes:  

![](/assets/brachy_smrna.png)

## Mapping

Map with bowtie, pipe to samtools to convert the output sam file to a bam file. Using the `-a` flag will tell bowtie to report all alignments.

```bash
bowtie -p 20 -q -n 0 -l 20 -a -S /home/tstuart/working_data_01/genomes/Bdistachyon/v2.1/assembly/BowtieIndex/bd21 adapertrim.fq \
| samtools view -b - > mapped.bam
```

## Post-processing

It's important to know the length of the smRNA, its mapping position in the genome, and whether it was uniquely mapped or not. First, we can convert the bam file to a bed file and add a column showing the smRNA length:

```bash
$ bedtools bamtobed -i mapped.bam \
| awk 'BEGIN {FS=OFS="\t"} {$7 = $3 - $2; print $0}' - | pigz - > mapped.bed.gz
```

Now we add a column showing the number of mapping positions for each smRNA (in R):

```r
library(dplyr)
library(data.table)
library(readr)

smrna <- fread("~/Downloads/mapped.bed", header = F, sep = "\t",
               col.names = c("chromosome", "start", "stop", "readName", "q", "strand", "length"))

smrna <- smrna %>%
  group_by(readName) %>%
  mutate(positions = n()) %>%
  ungroup() %>%
  select(-q) %>%
  mutate(unique = ifelse(positions > 1, FALSE, TRUE)) %>%
  arrange(chromosome, start)

write_tsv(smrna, "~/Desktop/brachy_smrna.tsv")
system("gzip ~/Desktop/brachy_smrna.tsv")
```

That's it! We now have a file that looks like this:

```bash
$ zcat brachy_smrna.tsv.gz | head
chromosome	start	stop	readName	strand	length	positions	unique
Bd1	3	27	SRR035616.1254270	-	24	143	FALSE
Bd1	5	26	SRR035616.1862694	-	21	153	FALSE
Bd1	6	21	SRR035616.3527272	-	15	218	FALSE
Bd1	10	34	SRR035616.1254270	-	24	143	FALSE
Bd1	12	33	SRR035616.1862694	-	21	153	FALSE
Bd1	13	28	SRR035616.3527272	-	15	218	FALSE
Bd1	17	41	SRR035616.1254270	-	24	143	FALSE
Bd1	19	40	SRR035616.1862694	-	21	153	FALSE
Bd1	20	35	SRR035616.3527272	-	15	218	FALSE
```
