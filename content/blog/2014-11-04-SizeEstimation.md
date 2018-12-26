---
title: Estimation of insert size from sam / bam files
date: 2014-11-04
---

I've been looking for a simple way to estimate paired-end insert sizes from mapped data. There are a few more complicated tools available (I think picard will do it), as well as simple scripts different people have written that don't exclude outliers. However outliers due to discordant read alignments will hugely change the mean, so need to be removed.

As I couldn't find anything that was a nice middle ground, I ended up writing my own simple python script that will take input from stdin and exclude outliers more that two standard deviations from the median, and print the mean and standard deviation to stdout.

It can be used to estimate the insert size from sam or bam files:

```
$ head -10000 mapped.sam | python mean_size.py
220 35
```

```
$ samtools view mapped.bam | head -10000 | python mean_size.py
220 35
```

<script src="https://gist.github.com/timoast/af73c0e9fac00187ee49.js"></script>
