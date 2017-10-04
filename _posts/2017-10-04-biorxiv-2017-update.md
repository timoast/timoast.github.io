---
title: bioRxiv 2017 update
author: Tim Stuart
date: '2017-10-04'
comments: true
layout: post
---

I first looked at the [biorxiv](https://www.biorxiv.org) submission data back in [March 2016](http://timoast.github.io/2016/03/01/biorxiv/). A lot has changed since then, and biorxiv has grown nearly 5-fold. Time for an update.

```r
collection_date <- ymd("2017_10_04")

dat <- fread("~/Documents/GitHub/biorxivData/data/biorxiv_data_2017_10_04.tsv") %>% 
  mutate(Age = collection_date - ymd(`Original submission`),
         Revised = `Original submission` != `Current submission`)
```

<!--break-->

### Submissions over time


```r
weekly <- dat %>%
  mutate(weeks_past = ceiling(Age / 7),
         `Submission week` = collection_date - weeks(weeks_past)) %>% 
  group_by(`Submission week`) %>%
  summarise(Submissions = n())

ggplot(weekly, aes(`Submission week`, Submissions)) +
  geom_point(stat = "identity") +
  geom_smooth() +
  ggtitle("bioRxiv submissions per week") +
  theme_bw()
```

![](/figure/2017-10-04-biorxiv-2017-update_files/figure-html/unnamed-chunk-3-1.png)<!-- -->

Last year the number of weekly submissions peaked at around 60, now it's 5x higher hitting 300 per week earlier in 2017.

How many of these submissions get revised?


```r
dat %>% 
  group_by(Revised) %>%
  summarise(n = n(), `%` = n/nrow(dat) * 100)
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["Revised"],"name":[1],"type":["lgl"],"align":["right"]},{"label":["n"],"name":[2],"type":["int"],"align":["right"]},{"label":["%"],"name":[3],"type":["dbl"],"align":["right"]}],"data":[{"1":"FALSE","2":"10603","3":"72.40014"},{"1":"TRUE","2":"4042","3":"27.59986"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

This is almost exactly the same percentage as I found last year.

### 2017 highlights

What have been the most popular preprints so far this year?


```r
days <- collection_date - ymd('2017-01-01')

dat %>%
  filter(Age < days) %>% 
  arrange(desc(`PDF views`)) %>% 
  head(10) %>% 
  select(Title)
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["Title"],"name":[1],"type":["chr"],"align":["left"]}],"data":[{"1":"Opportunities And Obstacles For Deep Learning In Biology And Medicine"},{"1":"Index Switching Causes “Spreading-Of-Signal” Among Multiplexed Samples In Illumina HiSeq 4000 DNA Sequencing"},{"1":"Regulation of Life Span by the Gut Microbiota in The Short-Lived African Turquoise Killifish"},{"1":"Sex Differences In The Adult Human Brain: Evidence From 5,216 UK Biobank Participants"},{"1":"The Reproducibility Of Research And The Misinterpretation Of P Values"},{"1":"Major flaws in \"Identification of individuals by trait prediction using whole-genome sequencing data\""},{"1":"The Beaker Phenomenon And The Genomic Transformation Of Northwest Europe"},{"1":"The Genomic History Of Southeastern Europe"},{"1":"The Human Cell Atlas"},{"1":"Comprehensive single cell transcriptional profiling of a multicellular organism by combinatorial indexing"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

Some really topical stuff (not surprisingly): p-value controversies, single cell genomics, the index-switching catastrophe, and the recent Venter debacle.

### Data

The data is available on my [github](https://github.com/timoast/biorxivData) to explore.

