---
title: "Fastest way to deduplicate pairs after reducing precision"
author: "Charles Plessy"
date: 'January 23, 2018'
output: 
  html_document: 
    keep_md: yes
---

To accelerate plotting and reduce the file size of vector graphics created with
the `plotCorrelation2` method, points are deduplicated after reducing their
precision to a smaller number of significant digits.  After benchmarking the
alternatives below, I chose the approach based on complex numbers, because their
performance is similar to `data.table` or `dplyr` objects, and because these
functions are available from standard R installations.




```r
set.seed(1)
pois1 <- rpois(1e7, 0.3)
pois2 <- rpois(1e7, 0.3)

nbin1 <- rnbinom(1e7, mu = 2, size = .1)
nbin2 <- rnbinom(1e7, mu = 2, size = .1)

# Same with a bit of noise
x     <- nbin1 * runif(1e7)
y     <- nbin2 * runif(1e7)

sx    <- signif(x, 2)
sy    <- signif(y, 2)
```

`DataFrame` tables are are coerced to `data.frame` objects by the
`duplicated.DataTable` function, therefore execution times are similar.


```r
f.df       <- function(x,y) {
  df <- unique(data.frame(x = x, y = y))
  data.frame(df$x, df$y)}

f.tbl <- function(x,y) {
  tbl <- dplyr::distinct(dplyr::tibble(x = x, y = y))
  data.frame(tbl$x, tbl$y)}

f.dt       <- function(x,y) {
  dt <- unique(data.table::data.table(x = x, y = y))
  data.frame(dt$x, dt$y)}

f.cplx     <- function(x,y) {
  u <- unique(complex(real=x, im=y))
  data.frame(Re(u), Im(u))}

f.cplx.Rle <- function(x,y) {
  u <- unique(S4Vectors::Rle(complex(real=x, im=y)))
  data.frame(Re(u), Im(u))}
```



```r
microbenchmark::microbenchmark(
  f.df       (pois1, pois2),
  f.dt       (pois1, pois2),
  f.tbl      (pois1, pois2),
  f.cplx     (pois1, pois2),
  f.cplx.Rle (pois1, pois2),  times = 10L)
```

```
## Unit: milliseconds
##                      expr       min        lq      mean    median        uq       max neval  cld
##        f.df(pois1, pois2) 5542.3731 5594.6264 5671.4823 5644.1065 5735.2541 5853.5424    10    d
##        f.dt(pois1, pois2)  285.1786  302.2330  329.5730  303.7755  308.4986  484.0237    10 a   
##       f.tbl(pois1, pois2)  247.2100  251.6843  279.2665  263.1273  275.9187  411.3219    10 a   
##      f.cplx(pois1, pois2)  668.0864  681.1732  733.4103  738.6093  786.9343  791.6873    10  b  
##  f.cplx.Rle(pois1, pois2)  798.5804  816.3618  924.7575  886.9679  913.6510 1459.9504    10   c
```


```r
microbenchmark::microbenchmark(
  f.df       (nbin1, nbin2),
  f.dt       (nbin1, nbin2),
  f.tbl      (nbin1, nbin2),
  f.cplx     (nbin1, nbin2),
  f.cplx.Rle (nbin1, nbin2),  times = 10L)
```

```
## Unit: milliseconds
##                      expr        min         lq       mean     median         uq        max neval   cld
##        f.df(nbin1, nbin2) 13724.8568 13795.2107 13844.4910 13825.8654 13880.0226 13996.2209    10     e
##        f.dt(nbin1, nbin2)  1431.5340  1438.4132  1476.0144  1448.3882  1463.9435  1607.0687    10    d 
##       f.tbl(nbin1, nbin2)   471.5175   473.9558   478.6101   475.1960   478.5863   508.0417    10 a    
##      f.cplx(nbin1, nbin2)   640.9068   652.0014   681.4527   657.1218   668.7779   797.4415    10  b   
##  f.cplx.Rle(nbin1, nbin2)   798.9246   804.1978   838.4483   816.2140   822.8612   956.9300    10   c
```


```r
microbenchmark::microbenchmark(
  f.df       (x, y),
  f.dt       (x, y),
  f.tbl      (x, y),
  f.cplx     (x, y),
  f.cplx.Rle (x, y),  times = 10L)
```

```
## Unit: milliseconds
##              expr        min         lq      mean    median        uq       max neval cld
##        f.df(x, y) 20419.0094 25717.6240 27444.943 26737.472 28482.954 37325.101    10   b
##        f.dt(x, y)  1943.3935  2014.8763  2472.712  2068.108  2171.388  4198.844    10  a 
##       f.tbl(x, y)  3600.8956  3633.6281  3867.547  3702.556  3791.214  4961.304    10  a 
##      f.cplx(x, y)   940.8565   997.6168  1576.777  1032.609  2792.510  2990.486    10  a 
##  f.cplx.Rle(x, y)  1085.7404  1097.5748  1373.175  1200.295  1710.728  2170.113    10  a
```


```r
microbenchmark::microbenchmark(
  f.df       (sx, sy),
  f.dt       (sx, sy),
  f.tbl      (sx, sy),
  f.cplx     (sx, sy),
  f.cplx.Rle (sx, sy),  times = 10L)
```

```
## Unit: milliseconds
##                expr        min         lq       mean     median         uq        max neval cld
##        f.df(sx, sy) 14004.1698 14085.9775 14411.9578 14376.2798 14665.1173 14964.3231    10   c
##        f.dt(sx, sy)  2112.3865  2153.9349  2212.3605  2212.3187  2254.7233  2324.2441    10  b 
##       f.tbl(sx, sy)   771.5590   795.3287   802.2588   807.4110   810.8238   825.1414    10 a  
##      f.cplx(sx, sy)   771.3897   777.3092   820.1100   806.8962   823.4534  1018.6156    10 a  
##  f.cplx.Rle(sx, sy)   887.5025   958.2897   957.9972   962.2164   971.3675  1024.1927    10 a
```
