---
title: "Fastest way to produce sorted abundances"
author: "Charles Plessy"
date: 'January 19, 2018'
output: 
  html_document: 
    keep_md: yes
---

Originally, the `.plotReverseCumulative` function was producing sorted abundances
using the data.table package, as in the function `f.dt` declared below.  After
benchmarking it with alternativesk, I switched to the approach implemented in
`f.Rle` because the data fed to `.plotReverseCumulative` is often already in
Rle format.




```r
set.seed(1)
pois <- rpois(1e7, 0.3)
nbin <- rnbinom(1e7, mu = 2, size = .1)
```


```r
f.df   <- function(x) as.data.frame(table(x))

f.dt   <- function(x) {
  v <- data.table::data.table(num = 1, nr_tags = x)
  v <- v[, sum(num), by = nr_tags]
  data.table::setkeyv(v, "nr_tags")
  as.data.frame(v)
}

f.Rle  <- function(x) {
  x <- S4Vectors::Rle(sort(x))
  data.frame(S4Vectors::runValue(x), S4Vectors::runLength(x))
}

f.rle  <- function(x) {
  x <- rle(sort(x))
  data.frame(x$values, x$lengths)
}

f.aggr <- function(x) {
  aggregate(x, by = list(x), FUN = length)
}
```


```r
microbenchmark::microbenchmark(f.df(pois), f.dt(pois), f.Rle(pois), f.rle(pois), f.aggr(pois), times = 10)
```

```
## Unit: milliseconds
##          expr       min        lq      mean    median        uq       max neval  cld
##    f.df(pois) 2547.3898 2576.0236 2632.0662 2649.5856 2673.0133 2703.1752    10   c 
##    f.dt(pois)  269.2576  277.5649  320.3839  294.6661  304.7185  521.1888    10 a   
##   f.Rle(pois)  250.1501  255.3788  330.7209  269.0381  282.4882  815.5093    10 a   
##   f.rle(pois)  516.7891  524.2573  560.9921  536.1125  612.9325  645.9420    10  b  
##  f.aggr(pois) 4180.2123 4290.8633 4382.2822 4373.8272 4433.9691 4628.4021    10    d
```


```r
microbenchmark::microbenchmark(f.df(nbin), f.dt(nbin), f.Rle(nbin), f.rle(nbin), f.aggr(nbin), times = 10)
```

```
## Unit: milliseconds
##          expr        min         lq       mean     median         uq        max neval  cld
##    f.df(nbin)  6677.9154  6736.9936  6826.6288  6776.7426  6920.7506  7131.7224    10   c 
##    f.dt(nbin)   851.7686   861.2842   905.3367   871.5668   988.5644   998.0199    10 a   
##   f.Rle(nbin)   847.3891   853.3732   864.6481   861.6981   872.6878   885.7835    10 a   
##   f.rle(nbin)  1147.6611  1152.8890  1197.0601  1162.5783  1263.9725  1295.2383    10  b  
##  f.aggr(nbin) 10407.4388 10584.1327 10695.4513 10647.4410 10763.4302 11188.3700    10    d
```
