---
title: "Fastest way to calculate BED12 blocks for UCSCData objects"
author: "Charles Plessy"
date: 'October 23, 2026'
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

Update (October 2026): building the blocks of all clusters at once
(`f.vectorised`) is about 700 × faster than looping on the clusters with
`f.mat.direct` (2 ms instead of 1.2 s for the 532 clusters of one sample), and
gives an identical result.  All the blocks are created with
a single call to `IRanges()` and split by cluster with `relist()` and
`PartitioningByEnd()`.  `exportToTrack` now uses this approach.




``` r
suppressPackageStartupMessages(library("GenomicRanges"))
suppressPackageStartupMessages(library("CAGEr"))
tc <- CAGEr::exampleCAGEexp |> CAGEr::tagClustersGR(1)
tcl <- split(tc, seq_along(tc))
decoded.mat <- rbind(width(tc), decode(tc$q_0.1), decode(tc$q_0.9))
```


``` r
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

# Vectorised over all clusters: one column per cluster in the matrices, up to
# three blocks per cluster (first base, interquantile range, last base), a single
# call to the IRanges constructor, and a split by cluster with relist().

f.vectorised <- function(decoded.mat) {
  width_value <- decoded.mat[1,]
  qLow_value  <- decoded.mat[2,]
  qUp_value   <- decoded.mat[3,]
  keep   <- rbind(qLow_value != 1, TRUE, qUp_value != width_value)
  starts <- rbind(1, qLow_value, width_value)[keep]
  ends   <- rbind(1, qUp_value,  width_value)[keep]
  relist(IRanges(starts, ends), PartitioningByEnd(cumsum(colSums(keep))))
}
```


``` r
tc[1]
```

```
## TagClusters object with 1 range and 7 metadata columns:
##     seqnames    ranges strand |           score    dominant_ctss   nr_ctss annotation  genes q_0.1 q_0.9
##        <Rle> <IRanges>  <Rle> |           <Rle>           <CTSS> <integer>      <Rle>  <Rle> <Rle> <Rle>
##   1    chr17  26050540      + | 22.310089406231 chr17:26050540:+         1   promoter grid1a     1     1
##   -------
##   seqinfo: 26 sequences (1 circular) from danRer7 genome
```

``` r
microbenchmark::microbenchmark(f.grline(tcl[[1]]), f.mat.direct(decoded.mat[,1]), f.mat.viaString(decoded.mat[,1]), times = 100)
```

```
## Unit: milliseconds
##                               expr      min       lq     mean   median       uq       max neval cld
##                 f.grline(tcl[[1]]) 6.651738 6.828249 7.343721 6.995603 7.414189 11.567481   100 a  
##     f.mat.direct(decoded.mat[, 1]) 2.133308 2.183652 2.331102 2.210591 2.266056  5.701192   100  b 
##  f.mat.viaString(decoded.mat[, 1]) 7.557194 7.739555 8.377135 7.897567 8.120161 19.510868   100   c
```

``` r
tc[5]
```

```
## TagClusters object with 1 range and 7 metadata columns:
##     seqnames            ranges strand |            score    dominant_ctss   nr_ctss annotation genes q_0.1 q_0.9
##        <Rle>         <IRanges>  <Rle> |            <Rle>           <CTSS> <integer>      <Rle> <Rle> <Rle> <Rle>
##   5    chr17 26453632-26453708      + | 1023.81202126363 chr17:26453667:+        16   promoter ttc7b    30    72
##   -------
##   seqinfo: 26 sequences (1 circular) from danRer7 genome
```

``` r
microbenchmark::microbenchmark(f.grline(tcl[[5]]), f.mat.direct(decoded.mat[,5]), f.mat.viaString(decoded.mat[,5]), times = 100)
```

```
## Unit: milliseconds
##                               expr      min       lq     mean   median       uq       max neval cld
##                 f.grline(tcl[[5]]) 7.491104 7.675015 8.007087 7.841324 8.004052 17.265202   100  a 
##     f.mat.direct(decoded.mat[, 5]) 2.836928 2.892983 2.987140 2.930143 2.983857  7.430684   100   b
##  f.mat.viaString(decoded.mat[, 5]) 7.543154 7.700290 7.901869 7.810298 7.959041 12.751089   100  a
```

All the clusters of one sample at once, which is what `exportToTrack` does:


``` r
length(tc)
```

```
## [1] 532
```

``` r
identical(IRangesList(apply(decoded.mat, 2, f.mat.direct)), f.vectorised(decoded.mat))
```

```
## [1] TRUE
```

``` r
microbenchmark::microbenchmark(
  loop       = IRangesList(apply(decoded.mat, 2, f.mat.direct)),
  vectorised = f.vectorised(decoded.mat),
  times = 10)
```

```
## Unit: milliseconds
##        expr         min         lq        mean      median          uq         max neval cld
##        loop 1231.833538 1234.47325 1285.285211 1238.099493 1239.874048 1717.153435    10  a 
##  vectorised    1.560775    1.66441    1.793331    1.775571    1.814277    2.419917    10   b
```
