---
title: "Fastest way to deduplicate pairs after reducing precision"
author: "Charles Plessy"
date: 'June 28, 2021'
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

`DataFrame` tables are coerced to `data.frame` objects by the
`duplicated.DataFrame` function, therefore execution times are similar.


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
##                      expr        min         lq      mean     median         uq        max neval
##        f.df(pois1, pois2) 6110.40949 6622.93829 9690.0623 8221.46085 8783.57322 26170.0297    10
##        f.dt(pois1, pois2)   60.02291   60.54951  133.4811   67.96144   74.58652   457.8002    10
##       f.tbl(pois1, pois2)  128.25698  132.04388  187.8383  133.13124  175.31787   524.5211    10
##      f.cplx(pois1, pois2)  405.69382  430.61868  676.9191  561.19835  783.51868  1332.7292    10
##  f.cplx.Rle(pois1, pois2)  505.85765  511.47287 1534.4948  731.14052 1255.44816  8163.3096    10
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
##                      expr       min        lq      mean    median        uq       max neval
##        f.df(nbin1, nbin2) 8052.2532 8441.2621 8807.6492 8772.0262 8970.0093 9651.1920    10
##        f.dt(nbin1, nbin2)  564.9376  590.3274  617.2036  602.0947  634.0411  713.5503    10
##       f.tbl(nbin1, nbin2)  228.5567  229.1026  252.5224  233.7307  241.0644  422.7830    10
##      f.cplx(nbin1, nbin2)  455.0756  460.2728  494.8075  487.4897  529.1486  540.6777    10
##  f.cplx.Rle(nbin1, nbin2)  550.5553  556.1877  716.5038  579.5182  667.8624 1247.6053    10
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
##              expr       min        lq      mean    median        uq        max neval
##        f.df(x, y) 8322.8257 8534.7328 9198.6651 9195.2933 9669.7323 10687.7408    10
##        f.dt(x, y)  726.0732  728.9718  843.5276  755.8585  785.2379  1473.6833    10
##       f.tbl(x, y)  582.6254  597.0590  610.0903  598.8600  605.4854   720.1217    10
##      f.cplx(x, y)  680.8757  686.3943  839.6481  701.8463  799.8486  1361.9492    10
##  f.cplx.Rle(x, y)  786.4886  795.6749  870.4042  805.1404  822.9168  1449.4373    10
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
##                expr       min        lq      mean    median        uq        max neval
##        f.df(sx, sy) 8155.3863 8260.2723 8828.2059 8759.4605 9053.9389 10242.4557    10
##        f.dt(sx, sy)  546.7252  580.4363  685.5099  596.9739  809.6086  1149.2348    10
##       f.tbl(sx, sy)  313.1412  317.8891  319.5706  319.3985  321.5431   326.8142    10
##      f.cplx(sx, sy)  573.0142  580.5656  614.5411  588.4361  638.7958   742.3095    10
##  f.cplx.Rle(sx, sy)  664.3078  673.3188  855.7151  744.8297  826.4204  1437.5688    10
```
