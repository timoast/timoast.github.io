---
title: Useful python code
layout: post
---
### Take command-line arguments

```python
from argparse import ArgumentParser

version = pkg_resources.require("program")[0].version
parser = ArgumentParser(description='program description')
group = parser.add_mutually_exclusive_group()
parser.add_argument('--version', action='version', version='%(prog)s '+str(version))
parser.add_argument('--option', help='option description', required=False, default=False, action='store_true')
parser.add_argument('-n', '--name', help='sample name', required=True)
```

## Lists

### List comprehension

If ever declaring and empty list then using `.append`, can be done with list comprehension.

```python
list3 = []
for x in list2:
    if x[0] in list1:
        list3.append(x)

list3 = [x for x in list2 if x[0] in list1]
```

Can also use this to read from a file in one line, eg to read column 4 in a file into a list in python:

```python
all_TEs = [line.rsplit()[4] for line in open("TAIR9_TE.bed", "r")]
```

Or read all lines into a list:

```python
all_lines = [line.strip("\n") for line in open("filename", "r")]
```

### Join items in list

```python
outfile.write('{l}\n'.format(l='\t'.join(list)))
```

If list items are not strings

```python
outfile.write('{l}\n'.format(l=map(str, list)))
```

Or using list comprehension

```python
stringVersion = [str(x) for x in inputList]
```

## Printing

### Print in colour

```python
class colour:
    PURPLE = '\033[95m'
    CYAN = '\033[96m'
    DARKCYAN = '\033[36m'
    BLUE = '\033[94m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    END = '\033[0m'

print colour.CYAN + 'This is blue' + colour.END
```

### Format text wrapping correctly in terminal

```python
import textwrap
print textwrap.dedent('print message')
```

## Strings

### Split strings by multiple delimiters

```python
accession_name = f.replace('calls_', '.').split('.')
```

## File reading / writing

### Read lines in tsv/csv file

```python
with open('merged_dmrs.bed', 'r') as f:
    for line in f:
        line = line.rsplit() # each item in line is now part of a list
```

## Find if coordinates overlap

```python
def overlap(start1, stop1, start2, stop2):
    """returns True if sets of coordinates overlap. Assumes coordinates are on same chromosome"""
    for y in xrange(start2, stop2+1):
        if start1 <= y <= stop1:
            return True
        else:
            pass
```

## Functions for mapping NGS reads

```python
import os
from subprocess import call


def peMap(proc):
    """
    paired-end reads
    """
    for filename in os.listdir('.'):
        if filename.endswith("_1.fastq"):
            sam = filename.split('_')
            sam = sam[0]
            call(["bowtie2", "-p{x}".format(x=proc), "-q", "-x /dd_stage/userdata/lister/data/genomes/bowtie2_indexes/tair9",
                  "-1 {f}".format(f=sam + '_1.fastq'), "-2 {f}".format(f=sam + '_2.fastq'), "-S {s}".format(s=sam + '.sam')])
        else:
            pass


def seMap(proc):
    """
    single-end reads
    """
    for filename in os.listdir('.'):
        if filename.endswith(".fastq"):
            sam = filename.split('.')
            sam = sam[0]
            call(["bowtie2", "-p{x}".format(x=proc), "-q", "-x /dd_stage/userdata/lister/data/genomes/bowtie2_indexes/tair9",
                  "-U {f}".format(f=filename), "-S {s}".format(s=sam + '.sam')])
        else:
            pass
```

## Split sra files into fastq

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
