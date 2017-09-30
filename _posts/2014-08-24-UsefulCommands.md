---
title: Useful Linux/Unix commands
layout: post
---
## go back to last directory

```
cd -
```

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

<!--break-->

## Split fasta file into a new file for each entry

```
awk '/>/{x="file_"++i;}{print > x".txt";}' input
```

## Split coordinate file into different files for each chromosome

```
awk  '{print > $1".txt"}' input
```

## Take command-line arguments

Required flags are followed by `:`

```bash
index=  proc=  path=  

while getopts x:pA opt; do
  case $opt in
  x)
      index=$OPTARG
      ;;
  p)
      proc=$OPTARG
      ;;
  A)
      path=${OPTARG%/}
      ;;
  esac
done
shift $((OPTIND - 1))
```

## Strip trailing slash from string

```bash
path=${OPTARG%/}
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

Add string prefix to a column, eg. add 'chr' to the first column

```bash
awk 'BEGIN {FS=OFS="\t"} {$1="chr"$1; print $0}' input_file.txt > output_file.txt
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

Compress recursively and store all files in a single compressed folder  

```bash
tar cvfz slam.tgz slam/
```

Compress recursively  

```bash
gzip -r directory/
```

Decompress files

```bash
tar -x -f file.tar
tar -x -z -f file.tgz
gunzip file.gz
```
