---
title: Useful Linux/Unix commands
layout: post
---
## Do something in a range
Move files

```bash
for i in $(seq 10); do mv chr$i ../genomes/; done
```

```bash
for directory in ./*; do
    if [ -d "$directory" ]; then
        cd $directory
        for myfile in $(ls -d *.sra);do
            mv $myfile /home/tstuart/working_data/1000genomes/$myfile
        done
        cd ..
    fi
done
```

## GNU Screen

Start screen: `screen -S [screen name]`  
List screens: `screen -list`  
Detach: `screen -d`  
Attach: `screen -r`  
Close screen: `ctr-a-d`  
Kill screen: `ctr-a :quit`  

## Human-readable path

```bash
echo -e ${PATH//:/'\n'}
```

## Searching

```bash
grep 'text' filename
grep 'text' file1 file2 file3
grep 'text1 text2' filename
grep --color 'text' filename
```

Search all files in a directory, show output in less  

```bash
grep -r 'text' /home/usr/ | less
```

Search for multiple strings  

```bash
egrep '(AT4G40030)|(AT4G40040)|(AT5G10980)' * > h3.3
```

## Pipe output into file
#### Overwrite contents of file

```bash
grep "text" filename.txt > output.txt
```

#### Append to file

```bash
grep "text" filename.txt >> output.txt
```

## Counting
Count lines in chr1 file

```bash
wc -l chr1
```

Count characters in chr1 file

```bash
wc -c chr1
```

## List

List and sort  

```bash
ls -ls
```

List and sort by size  

```bash
ls -ls -S
```

Count number of directories in current directory  

```bash
ls -l | grep ^d | wc -l
```

## Sorting

Sort file by numerical order of first column, save as sorted_list.txt  

```bash
sort -nk 1 list.txt > sorted_list.txt
```

Sort file by order of first column, then numerical order of second column, save as sorted_list.txt  

```bash
sort -k1,1 -nk2,2 list.txt > sorted_list.txt
```

Sort descending  

```bash
sort -rn -k3 file.txt > sorted.txt
```

overwrite  

```bash
sort -rn -k3 -o file.txt file.txt
```

## Rows
Delete first line

```bash
sed 1d file.txt > headerless_file.txt
```

Delete lines 1-3 inclusive  

```bash
sed 1,3d file.txt > newfile.txt
```

Delete lines containing string  

```bash
sed '/ATC/d' file.txt > file_mod.txt
```

Overwrite file  

```bash
sed -i.bak 1d file.txt
```

## Columns
Write columns to new file. Can also reorder columns.  

```bash
awk 'BEGIN {FS=OFS="\t"} {print $2,$4,$5,$6,$7,$8,$1,$3}' file.txt > outfile.txt
awk 'BEGIN {FS=OFS="\t"} {print $7,$1,$8,$2,$3,$4,$5,$6}' cmt2_targets.tsv > cmt2_targets_ordered.tsv
```

Add column  

```bash
awk '{print $0, "cmt2"}' cmt2_targets_ordered.tsv > newfile.tsv
```

Remove first column from all files  

```bash
for myfile in $(ls);do
    awk 'BEGIN {FS=OFS="\t"} {$1="";sub("\t","")}1' $myfile > mod_$myfile
done
```

Search only in one column

```bash
awk '{if ($1 == 1) print $2}' p1 > chr1_p1
```

## Joining files

Join files with each starting on a new line  

```bash
for filename in *.fa; do
    cat "${filename}"
    echo
done > output.fa
```
Merge files of the same format

```bash
cat p1_* | sort -nk1,1 -k2,2 | uniq > p1
```

## Comparing Files  

Files should first be sorted


Find lines that are common or different between two files  

```bash
comm -i file1 file2 > output.txt
```

* `-i` is case-insensitive  
* outputs three columns. First is lines only in file1, second is lines only in file2, third is lines common in both.

## Compressing files  

Compress recursively  

```bash
tar cvfz slam.tgz slam/
```

Decompress files

```bash
tar -x -f file.tar
tar -x -z -f file.tgz
```

## Download data from SRA  

### Using wget  

```
wget -r --no-parent --reject "index.html*" ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP/SRP005/SRP005399/
```
Alternatively, you can download a list of accessions from the SRA website and use that file to call `wget` for all the accessions.  

```python
from subprocess import call


with open('SraAccList.txt', 'r') as accessions:
    for row in accessions:
        acc = row.strip('\n')
        print "Downloading {a}".format(a=acc)
        call(["wget",
              "-r",
              "--no-parent",
              "ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR492/{acc}/*".format(acc=acc)])
```

#### Split into fastq files

```bash
fastq-dump --split-3 ./SRP005399/*
```

`--split-3` causes PE files to be split into `*_1.fastq` and `*_2.fastq`

bash loop  

```bash
for directory in ./*; do
    if [ -d "$directory" ]; then
        cd $directory
        for myfile in $(ls -d *.sra);do
            fastq-dump --split-3 -v $myfile
        done
        cd ..
    fi
done
```

Python function  

```python
import os
from subprocess import call


def fastqSplit():
    for filename in os.listdir('.'):
        if filename.endswith('.sra'):
            print 'processing {n}'.format(n=filename)
            call(['fastq-dump', '--split-3', '-v', filename])
        else:
            pass
```

### Using SRA toolkit  

SRA toolkit doesn't seem to work well for downloading bulk data. `wget` is a much better option as it allows download of whole folders, gives more descriptive output as it goes. SRA often fails and sometimes gives no error. Download with `wget` and use `fastq-split` to get fastq files.  

```bash
fastq-dump SRR534224 &
```

## RNA-seq

### Mapping  

```
tophat -p8 -G /home/lister/working_data/data/genomes/annotations/tair10/TAIR10_gene_TE_illumina.gtf --transcriptome-index=transcript_index/tair10 /home/lister/working_data/data/genomes/bowtie2_indexes/tair9 SRR501604.fastq,SRR501605.fastq
```
* `--transcriptome-index=transcript_index/tair10` creates transcriptome index file saved in transcript_index/ with name `tair10.*`. This can be reused if mapping to the same transcriptome (ie. gff file)  
* `-p8` uses 8 cores  
* `-G` specifies gff file. Optional  
* If using PE reads check manual  

## Bedtools  

[Good tutorial from the Quinlan lab](http://quinlanlab.org/tutorials/cshl2013/bedtools.html)  

Need to order files so that chromosome, start, stop and first three columns (BED format). Also need to remove header. Can do these steps with `awk` and `sed 1d`.  

### Merge  

```bash
bedtools merge -i dmrs_ddc.tsv -d 100 > dmrs_ddc_merged.tsv
```
* `-d` is distance between coordinates that will still be merged  
* `-i` is input file  

### Intersect  

```bash
bedtools intersect -a tair10_tes.txt -b dmrs_ddc_merged.tsv -wa -f 0.50 > ddc_targets.tsv
```
* `-a` is file a  
* `-b` is file b  
* `-wa` is write a  
* `-f` is minimum overlap percentage  

## Bowtie  

### Mapping PE data  

```
bowtie2 -p8 --local --fr -q -R 5 -N 1 -x [path/to/bowtie2/index] -X 1000 -1 [mate1.fq] -2 [mate2.fq] -S [alignment.sam] --un-conc ./discordant/
```
* `-p` is number of processors  
* `--local` allows soft-clipping of reads to improve alignment  
* `-q` specifies that reads are in fastq format  
* `-x` gives path to index  
* `-S` gives name of output file  
* `-R` is number of times bowtie will try to 're-seed' repetitive seeds. Default 2  
* `-N` is number of mismatches allowed in seed. 0 or 1.  
* `--no-mixed` tells bowtie to find alignments only when both pairs can be aligned.  
* `--fr` means mate pairs are ordered in forward then reverse orientation. Can do `--ff`, `--rf`.  
* `-X` specifies maximum insert size. Default 500.  
* `-I` specifies minimum insert size. Default 0 (no minimum).  
* `--un-conc` specifies path and file to write discordant alignments to. Note that these are just the `fastq` reads, not alignments.

## Samtools

samtools can read from a stream, so can pipe output in from other tools (eg bowtie to get `.bam` output) or other samtools commands.

### Convert from `sam` to `bam`
```
samtools view -bS file.sam > file.bam
```

### Sort bamfile
```
samtools sort file.bam sorted
```
* outputs `sorted.bam` file

### Samblaster

Extracts reads from `.sam` alignment files

#### Get discordant reads
```
bowtie2 [options] | samblaster -e -d file.sam | samtools view -bS - > file.bam
```
* Print output from bowtie to stderr to pipe directly into samblaster
* Saves having to search through sam or bam file later to extract discordant reads
* Pipe output directly into `samtools view` to save as `.bam` file (these won't be the reads in the samblaster output)
