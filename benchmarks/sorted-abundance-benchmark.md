---
title: "Fastest way to produce sorted abundances"
author: "Charles Plessy"
date: 'June 28, 2021'
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
##          expr       min        lq      mean    median        uq       max neval
##    f.df(pois)  415.8674  439.0261  475.4740  447.0769  508.7650  608.3248    10
##    f.dt(pois)  186.5875  197.6726  240.8832  242.7743  260.1480  368.5690    10
##   f.Rle(pois)  101.2502  102.0797  160.9364  105.4015  127.5486  571.2380    10
##   f.rle(pois)  197.3088  210.7383  221.8933  220.4917  225.4271  264.0947    10
##  f.aggr(pois) 2125.0768 2191.2064 2259.0070 2247.5334 2330.3711 2420.9537    10
```


```r
microbenchmark::microbenchmark(f.df(nbin), f.dt(nbin), f.Rle(nbin), f.rle(nbin), f.aggr(nbin), times = 10)
```

```
## Unit: milliseconds
##          expr       min        lq      mean    median        uq       max neval
##    f.df(nbin) 4545.5276 4642.2336 4674.8707 4661.2145 4692.9570 4857.7300    10
##    f.dt(nbin)  335.1644  369.9447  405.9829  382.1086  409.7706  608.3491    10
##   f.Rle(nbin)  424.9090  456.2392  470.5724  476.7369  490.5005  500.1160    10
##   f.rle(nbin)  577.1908  591.7957  623.9684  610.1879  623.4713  720.8229    10
##  f.aggr(nbin) 6632.6367 6716.8786 6836.4560 6843.9143 6867.1097 7230.9565    10
```
