---
title: Useful bioinformatics
author: Tim Stuart
date: '2017-09-14'
comments: true
layout: post
---

## Trim reads

With [cutadapt](http://cutadapt.readthedocs.io/en/stable/)

```bash
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -m 30 -o output.fq.gz input.fq.gz
```

With [trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)

```bash
trimmomatic PE -threads 10 \
      reads_1.fastq.gz reads_2.fastq.gz \
      reads_1_trim.fq reads_1_se_trim.fq reads_2_trim.fq reads_2_se_trim.fq \
      ILLUMINACLIP:TruSeq2-PE.fa:2:30:10 \
      LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```

<!--break-->

## Filter short reads from FASTQ

```bash
gzip -d -c *.fastq.gz | paste - - - - \
  | awk 'length($2) > 25' \
  | tr "\t" "\n" 
```

## Map RNA-seq

```bash
# build index
STAR --runThreadN 10 \
     --runMode genomeGenerate \
     --genomeDir genome \
     --genomeFastaFiles genome.fa \
     --sjdbGTFfile genes.gtf
     
# map
STAR --runThreadN 20 \
     --genomeDir genome \
     --alignIntronMax 5000 \
     --alignIntronMin 10 \
     --readFilesCommand zcat \
     --quantMode GeneCounts \
     --outSAMtype BAM SortedByCoordinate \
     --outFileNamePrefix output_ \
     --readFilesIn input.fq.gz
```

## Pseudomap RNA-seq with kallisto

```bash
# build index 
kallisto index -i index.idx Arabidopsis_thaliana.TAIR10.cds.all.fa

#quantify
for filename in *.fq.gz; do
    samplename=(${filename//.fq.gz/ })
    kallisto quant -i index.idx -o $samplename -b 100 --single -l 150 -s 20 $filename -t 4
done
```

## Create coverage track from bam file

`bamCoverage` is part of [deepTools](http://deeptools.readthedocs.io/en/latest/index.html).

```bash
bamCoverage -b reads.bam -o coverage.bw -p 10
```

## Create UCSC browser track

You need a track hub, a trackdb file, and track files. These should be placed in a publicly accessible place.

Track hub:

```
hub Arabidopsis
shortLabel RootRNAseq
longLabel Root
genomesFile genomes.txt
email timstuart90@gmail.com
```

trackdb

```
include athTha1.annotationTracks.txt
include root/root.track.txt
```

Track

```
track root_coverage
type bigWig
bigDataUrl protoplast/root.bw
shortLabel root
longLabel root_coverage
visibility squish
priority 2
maxHeightPixels 100:60:10
```

On UCSC go to my data > track hubs > my hubs and enter the address of the track hub file.

## Bedtools  

```bash
bedtools merge -i dmrs_ddc.tsv -d 100 > dmrs_ddc_merged.tsv
```
* `-d` is distance between coordinates that will still be merged  
* `-i` is input file  

```bash
bedtools intersect -a a.bed -b b.bed -wa -f 0.50 > targets.bed
```

* `-wa` writes `-a` coordinates  
* `-f` is minimum overlap percentage  

## Bowtie  

### Mapping PE data  

```
bowtie2 -p8 --local --fr -q -R 5 -N 1 -x [path/to/bowtie2/index] -X 1000 \
        -1 read_1.fq.gz -2 read_2.fq.gz -S aligned.sam 
```
* `-p` is number of processors  
* `--local` allows soft-clipping of reads to improve alignment  
* `-q` specifies that reads are in fastq format  
* `-x` gives path to index  
* `-S` gives name of output file  
* `-R` is number of times bowtie will try to 're-seed' repetitive seeds. Default 2  
* `-N` is number of mismatches allowed in seed. 0 or 1.  
* `--fr` means mate pairs are ordered in forward then reverse orientation. Can do `--ff`, `--rf`.  
* `-X` specifies maximum insert size. Default 500.  

## Samtools

samtools can read from a stream, so can pipe output in from other tools (eg bowtie to get `.bam` output) or other samtools commands.

### Convert from `sam` to `bam`
```
samtools view -bS file.sam > file.bam
```

### Sort bamfile

```
samtools sort -@ 20 -O bam -T tmp file.bam 
```

### Samblaster

Extracts reads from `.sam` alignment files

#### Get discordant reads
```
bowtie2 [options] | samblaster -e -d file.sam | samtools view -bS - > file.bam
```
* Print output from bowtie to stderr to pipe directly into samblaster
* Saves having to search through sam or bam file later to extract discordant reads
* Pipe output directly into `samtools view` to save as `.bam` file (these won't be the reads in the samblaster output)
