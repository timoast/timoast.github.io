---
title: Extract reads from bam file by read name
layout: post
---

While there are very fast and easy ways to extract reads from a bam file according to mapping location, extracting reads by read name is more difficult.

Simple methods, like using `grep`, are incredibly slow if you want to look for more than a few reads.

Luckily, `pysam` allows you to index a bam file by read name (using `pysam.IndexedReads(AlignmentFile)`) while keeping the sort order. However, the index is stored in memory so this can use a lot of RAM (my tests with a 5.7 GB bam file used about 9 GB RAM).

I wrote a small python script (below) that uses pysam to extract reads by read name from a bam file.

Extracting 10 reads from a 5.7 GB bam file, just using `grep` is slightly faster than the python script:

```shell
timstuart Altai-5$  time samtools view Altai-5_filtered.bam | grep -f reads.txt > extracted

real    2m10.088s
user    2m23.107s
sys 0m35.470s

timstuart Altai-5$  time python extract_reads.py -b Altai-5_filtered.bam -n reads.txt -o python_extracted.bam

real    3m36.444s
user    3m14.744s
sys 0m18.867s

```

But extracting 200 reads using python is fast, while using `grep` ran for over an hour before I stopped the process.

```shell
timstuart Altai-5$  time python extract_reads.py -b Altai-5_filtered.bam -n reads.txt -o python_extracted.bam

real    3m30.015s
user    3m12.738s
sys 0m15.055s
```


<script src="https://gist.github.com/timoast/2264a79f93b3f1cb3aac.js"></script>
