---
title: "Fastest way to find dominant CTSS"
author: "Charles Plessy"
date: 'February 21, 2022'
output: 
  html_document: 
    keep_md: yes
editor_options: 
  chunk_output_type: console
---

# Purpose

We use `data.table` to find the _dominant CTSS_ in _tag clusters_.  I am
implementing a similar function for _consensus clusters_.  Shall I do it
with Bioconductor core functions, `data.table` or `dplyr`?

# Setup

We have CTSSes and clusters.  The clusters may overlap with each other.


```r
library(CAGEr) |> suppressPackageStartupMessages()
ce <- exampleCAGEexp
ctss <- CTSScoordinatesGR(ce)
score(ctss) <- CTSSnormalizedTpmDF(ce) |> DelayedArray::DelayedArray() |> rowSums()
clusters <- consensusClustersGR(ce)
clusters$annotation <- clusters$genes <- clusters$exprClass <- NULL # Simplify
```

The functions are given a lookup table associating CTSSes with clusters, and
a function to find the dominant CTSS and break ties by selecting the one
at the center.


```r
o <- findOverlaps(clusters, ctss)

find.dominant.idx <- function (x) {
  # which.max is breaking ties by taking the last, but this will give slightly
  # different biases on plus an minus strands.
  w <- which(x == max(x))
  w[ceiling(length(w)/2)]
}
```

## _Bioconductor_ way


```r
bioC <- function () {
  s <- extractList(score(ctss), o)
  m <- sapply(s, find.dominant.idx)
  grl <- extractList(granges(ctss), o)
  dom <- mapply(`[`, grl, m) |> GRangesList() |> unlist()
  clusters$dominant_ctss <- dom
  clusters$tpm.dominant_ctss <- max(s)
  clusters
}

bioC() -> bioC_results
```

## _Bioconductor_ plus base `R`


```r
bioC2 <- function () {
  rl <- rle(queryHits(o))$length
  cluster_start_idx <- cumsum(c(1, head(rl, -1))) # Where each run starts
  grouped_scores <- extractList(score(ctss), o)
  local_max_idx <- sapply(grouped_scores, find.dominant.idx) -1  # Start at zero
  global_max_ids <- cluster_start_idx + local_max_idx
  clusters$dominant_ctss <- granges(ctss)[subjectHits(o)][global_max_ids]
  clusters$tpm.dominant_ctss <- score(ctss)[subjectHits(o)][global_max_ids]
  clusters
}

bioC2() -> bioC_results2
```

## `data.table` way


```r
library("data.table") |> suppressPackageStartupMessages()
dataTable <- function() {
  dt <- ctss |> as.data.frame() |> data.table::as.data.table()
  dt$id <- dt$cluster |> as.factor() |> as.integer()
  dom <- dt[ , list( seqnames[1]
                   , strand[1]
                   , pos[find.dominant.idx(score)]
                   , sum(score)
                   , max(score))
            , by = id]
  setnames(dom, c( "cluster", "chr", "strand"
                 , "dominant_ctss", "tpm", "tpm.dominant_ctss"))
  
  # Let's be really sure we do not mix up things
  stopifnot(identical(clusters$consensus.cluster, dom$cluster))
  clusters$dominant_ctss <- GRanges(dom$chr, dom$dominant_ctss, dom$strand)
  seqinfo(clusters$dominant_ctss) <- seqinfo(ctss)
  clusters$tpm.dominant_ctss <- dom$tpm.dominant_ctss
  clusters  
}

dataTable() -> dataTable_results
```

## `dplyr` way


```r
dplyr <- function() {
  tb <- ctss |> as.data.frame() |> tibble::as_tibble()
  tb$id <- tb$cluster |> as.factor() |> as.integer()
  dom <- tb |> dplyr::group_by(id) |>
      dplyr::summarise(seqnames = unique(seqnames), strand = unique(strand),
                       tpm.dominant_ctss = max(score), dominant_ctss = pos[find.dominant.idx(score)])
  clusters$dominant_ctss <- GRanges(dom$seqnames, dom$dominant_ctss, dom$strand)
  seqinfo(clusters$dominant_ctss) <- seqinfo(ctss)
  clusters$tpm.dominant_ctss <- dom$tpm.dominant_ctss
  clusters
}

dplyr() -> dplyr_results
```

# Checks and benchmark


```r
bioC_results
```

```
## ConsensusClusters object with 805 ranges and 5 metadata columns:
##                             seqnames            ranges strand |       score
##                                <Rle>         <IRanges>  <Rle> |   <numeric>
##            chr17:26027430:+    chr17          26027430      + |     64.1719
##            chr17:26050540:+    chr17          26050540      + |     22.3101
##            chr17:26118088:+    chr17          26118088      + |     64.1719
##            chr17:26142853:+    chr17          26142853      + |     56.0704
##            chr17:26166954:+    chr17          26166954      + |     64.1719
##                         ...      ...               ...    ... .         ...
##   chr17:32706021-32706407:+    chr17 32706021-32706407      + | 136242.4564
##            chr17:32706605:+    chr17          32706605      + |     64.1719
##   chr17:32707132-32707170:+    chr17 32707132-32707170      + |    245.4415
##   chr17:32707322-32707376:+    chr17 32707322-32707376      + |    357.4316
##   chr17:32708847-32708958:+    chr17 32708847-32708958      + |   8395.7855
##                             consensus.cluster         tpm    dominant_ctss
##                                     <integer>   <numeric>        <GRanges>
##            chr17:26027430:+                 1     64.1719 chr17:26027430:+
##            chr17:26050540:+                 2     22.3101 chr17:26050540:+
##            chr17:26118088:+                 3     64.1719 chr17:26118088:+
##            chr17:26142853:+                 4     56.0704 chr17:26142853:+
##            chr17:26166954:+                 5     64.1719 chr17:26166954:+
##                         ...               ...         ...              ...
##   chr17:32706021-32706407:+               801 136242.4564 chr17:32706231:+
##            chr17:32706605:+               802     64.1719 chr17:32706605:+
##   chr17:32707132-32707170:+               803    245.4415 chr17:32707135:+
##   chr17:32707322-32707376:+               804    357.4316 chr17:32707322:+
##   chr17:32708847-32708958:+               805   8395.7855 chr17:32708890:+
##                             tpm.dominant_ctss
##                                     <numeric>
##            chr17:26027430:+           64.1719
##            chr17:26050540:+           22.3101
##            chr17:26118088:+           64.1719
##            chr17:26142853:+           56.0704
##            chr17:26166954:+           64.1719
##                         ...               ...
##   chr17:32706021-32706407:+        14993.7493
##            chr17:32706605:+           64.1719
##   chr17:32707132-32707170:+           99.6448
##   chr17:32707322-32707376:+           67.2963
##   chr17:32708847-32708958:+         1319.9060
##   -------
##   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```

```r
identical(bioC_results, bioC_results2)
```

```
## [1] TRUE
```

```r
identical(bioC_results, dataTable_results)
```

```
## [1] TRUE
```

```r
identical(bioC_results, dplyr_results)
```

```
## [1] TRUE
```

```r
microbenchmark::microbenchmark(bioC(), bioC2(), dataTable(), dplyr())
```

```
## Unit: milliseconds
##         expr        min         lq       mean     median         uq        max
##       bioC() 3455.15544 3588.41609 3730.71735 3660.22403 3832.07402 4245.39739
##      bioC2()   34.00017   36.20774   38.40580   37.57237   39.03414   50.95568
##  dataTable()   58.93020   62.93159   67.10066   64.74682   69.38356   90.61266
##      dplyr()  107.36517  112.24707  119.96521  115.04563  123.53149  185.93090
##  neval cld
##    100   c
##    100 a  
##    100 a  
##    100  b
```

The approach combining `findOverlaps` from _Bionductor_ and cumulative sums from
base _R_ works the best on test data, with comparable performance with `dplyr`
and `data.table`.  This opens the way to the removal of the dependency to
`data.table`.  In the future, if we start to import `dplyr` for reasons related
to `ggplot2` or `plyranges`, we may might replace the _Bioc + base R_ version
with a `dplyr` one, that is probably easier to read for most contributors. 
