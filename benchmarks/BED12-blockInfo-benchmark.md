---
title: "Fastest way to calculate BED12 blocks for UCSCData objects"
author: "Charles Plessy"
date: 'June 29, 2021'
output: 
  html_document: 
    keep_md: yes
editor_options: 
  chunk_output_type: console
---

In the [BED12 format](https://genome.ucsc.edu/FAQ/FAQformat.html#format1),
the position of each box in the object is encoded in the `blockCount`,
`blockSize` and `blockStarts` columns.

In `UCSCData` track objects of the _rtracklayer_ package, the information is
encoded as an `IRangesList` metadata column, where each element has one
`IRanges` object representing the blocks.

In _CAGEr_'s `TagClusters` objects, the position of the central blocks is
represented by the quantile information columns.

The problem here is to create this `IRanges` list in a computationally efficient
way in _CAGEr_'s `exportToTrack` function.


The best approach, `f.mat.direct` is not exactly fast, but is 5 × faster than my
original attempt (`f.grline`).




```r
suppressPackageStartupMessages(library("GenomicRanges"))
suppressPackageStartupMessages(library("CAGEr"))
tc <- CAGEr::exampleCAGEexp |> CAGEr::tagClustersGR(1)
tcl <- split(tc, seq_along(tc))
decoded.mat <- rbind(width(tc), decode(tc$q_0.1), decode(tc$q_0.9))
```


```r
# Operate directly on each line

f.grline <- function(grline) {
  # A GRanges line is a GRanges object of length 1
  # This function computes a blocks value for each line
  qLow_value <- decode(grline$q_0.1)
  qUp_value  <- decode(grline$q_0.9)
  ir <- IRanges()                      |>
    c( if(qLow_value != 1) IRanges(1)) |>
    c( IRanges(qLow_value, qUp_value)) |>
    c( if(qUp_value != width(grline)) IRanges(width(grline)))
  ir
}

# Decode the Rle objects once for all and loop on each entry

f.mat.direct <- function(x) {
  width_value <- x[1]
  qLow_value <- x[2]
  qUp_value  <- x[3]
  ir <- IRanges()                      |>
    c( if(qLow_value != 1) IRanges(1)) |>
    c( IRanges(qLow_value, qUp_value)) |>
    c( if(qUp_value != width_value) IRanges(width_value))
}

# Same but try to call the IRanges constructor only once.

f.mat.viaString <- function(x) {
  width_value <- x[1]
  qLow_value <- x[2]
  qUp_value  <- x[3]
  str <- c(
    c( if(qLow_value != 1) "1"),
    c( paste0(qLow_value, "-", qUp_value)),
    c( if(qUp_value != width_value) width_value)
  )
  IRanges(str)
}
```


```r
tc[1]
```

```
## TagClusters object with 1 range and 6 metadata columns:
##     seqnames    ranges strand |           score   nr_ctss dominant_ctss tpm.dominant_ctss q_0.1 q_0.9
##        <Rle> <IRanges>  <Rle> |           <Rle> <integer>     <integer>             <Rle> <Rle> <Rle>
##   1    chr17  26050540      + | 22.310089406231         1      26050540   22.310089406231     1     1
##   -------
##   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```

```r
microbenchmark::microbenchmark(f.grline(tcl[[1]]), f.mat.direct(decoded.mat[,1]), f.mat.viaString(decoded.mat[,1]), times = 100)
```

```
## Unit: milliseconds
##                               expr      min        lq      mean    median        uq       max neval
##                 f.grline(tcl[[1]]) 9.178579  9.645029 10.704038 10.284728 11.355642 15.512259   100
##     f.mat.direct(decoded.mat[, 1]) 2.810769  2.961396  3.528324  3.302803  3.606683  9.024574   100
##  f.mat.viaString(decoded.mat[, 1]) 9.342554 10.126850 11.103980 10.933802 11.565727 22.582231   100
```

```r
tc[5]
```

```
## TagClusters object with 1 range and 6 metadata columns:
##     seqnames            ranges strand |            score   nr_ctss dominant_ctss tpm.dominant_ctss q_0.1 q_0.9
##        <Rle>         <IRanges>  <Rle> |            <Rle> <integer>     <integer>             <Rle> <Rle> <Rle>
##   5    chr17 26453632-26453708      + | 1023.81202126363        16      26453667  288.917306039892    30    72
##   -------
##   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```

```r
microbenchmark::microbenchmark(f.grline(tcl[[5]]), f.mat.direct(decoded.mat[,5]), f.mat.viaString(decoded.mat[,5]), times = 100)
```

```
## Unit: milliseconds
##                               expr       min        lq      mean    median        uq       max neval
##                 f.grline(tcl[[5]]) 10.015783 10.576726 11.313155 11.094960 12.048345 14.515133   100
##     f.mat.direct(decoded.mat[, 5])  3.665505  3.838581  4.120408  3.974128  4.377097  5.110636   100
##  f.mat.viaString(decoded.mat[, 5])  9.529413 10.062541 10.744386 10.541201 11.260659 13.957480   100
```
