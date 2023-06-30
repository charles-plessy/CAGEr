---
title: "rowSums on Rle DataFrames"
author: "Charles Plessy"
date: "2023-06-30"
output: 
  html_document: 
    keep_md: yes
editor_options: 
  chunk_output_type: console
---

## Purpose

There is no native `rowSums` function for `DataFrame`s of `Rle`-encoded
numerical values.  What is the best strategy for an implementation ?

## Setup


```r
library("CAGEr")   |> suppressPackageStartupMessages()
library("ggplot2") |> suppressPackageStartupMessages()

# DF <- CTSStagCountDF(exampleCAGEexp)
DF <- ZebrafishDevelopmentalCAGE::ZebrafishCAGE |> assay()

f_Reduce_native         <- function(DF) Reduce(`+`, DF, Rle(0))
f_decode_reduce         <- function(DF) Rle(Reduce(`+`, lapply(DF, decode), 0))
f_rowSums_as.data.frame <- function(DF) Rle(rowSums(as.data.frame(DF)))
f_rowSums_decode        <- function(DF) Rle(rowSums(data.frame(lapply(DF, decode))))
f_rowSums_delayedarray  <- function(DF) Rle(rowSums(DelayedArray::DelayedArray(DF)))
```

## Benchmark


```r
(microbench_out <- microbenchmark::microbenchmark(times = 100,
  f_Reduce_native(DF),
  f_decode_reduce(DF),
  f_rowSums_as.data.frame(DF),
  f_rowSums_decode(DF),
  f_rowSums_delayedarray(DF)))
```

```
## Unit: milliseconds
##                         expr       min        lq      mean    median        uq
##          f_Reduce_native(DF) 3123.5140 3320.7487 3532.0372 3451.9774 3615.9572
##          f_decode_reduce(DF)  311.3833  355.1572  494.7205  387.1896  452.8512
##  f_rowSums_as.data.frame(DF)  448.3442  520.1004  671.9624  572.1575  652.3051
##         f_rowSums_decode(DF)  438.3752  505.2630  618.4694  559.2740  622.9931
##   f_rowSums_delayedarray(DF) 2819.1220 3206.4391 3727.6560 3750.7067 4100.5520
##       max neval  cld
##  4777.751   100   c 
##  1516.292   100 a   
##  1567.608   100  b  
##  1624.716   100 ab  
##  5138.834   100    d
```

```r
# https://statisticsglobe.com/microbenchmark-package-r
ggplot(microbench_out, aes(x = time, y = expr, color = expr)) +  # Plot performance comparison
  geom_boxplot() + 
  scale_x_log10("time (miliseconds)")
```

![](rowSums_on_Rle_DF_files/figure-html/benchmark-1.png)<!-- -->

## Result

The winner is:

```
f_decode_reduce <- function(DF) Rle(Reduce(`+`, lapply(DF, decode), 0))
```

Explained with more verbose code:

```
f_decode_reduce <- function(DF) {
  list_of_numerical_vectors <- lapply(DF, decode)
  parallel_sum              <- Reduce(`+`, list_of_numerical_vectors, 0)
  result                    <- Rle(parallel_sum)
  return(result)
}
```
