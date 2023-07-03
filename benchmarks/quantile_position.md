---
title: "Quantile positions in cumulative sums"
author: "Charles Plessy"
date: "2023-06-30"
output: 
  html_document: 
    keep_md: yes
editor_options: 
  chunk_output_type: console
---

## Purpose

The function to determine quantile positions is one of the slowest in CAGEr.
Can we speed it up?

## Setup


```r
library("CAGEr")   |> suppressPackageStartupMessages()
library("ggplot2") |> suppressPackageStartupMessages()
```

### Example data


```r
cum.sums <- CTSScumulativesTagClusters(exampleCAGEexp)[[1]]
clusters <- tagClustersGR(exampleCAGEexp)[[1]]
```

## Current implementation at the time the benchmark was written.


```r
getQuantilepos_Vectorized_tail_decoded <- Vectorize(vectorize.args = "cum.sum", function(q, cum.sum) {
  cum.sum <- decode(cum.sum) # A microbenchmark showed it it 3 times faster when applying decode() now
  c.max <- tail(cum.sum,1) # Max is last element since x is a cumulative sums.
  treshold <- c.max * q
  which.max(cum.sum >= treshold)
})
getQuantilepos_Vectorized_tail_decoded(.9, cum.sums) |> head()
```

```
##  1  2  3  4  5  6 
##  1  1  1  1 72  1
```

## Alternative functions.


```r
getQuantilepos_Vectorized_tail <- Vectorize(vectorize.args = "cum.sum", function(q, cum.sum) {
  c.max <- tail(cum.sum,1) # Max is last element since x is a cumulative sums.
  treshold <- c.max * q
  which.max(cum.sum >= treshold)
})
getQuantilepos_Vectorized_tail(.9, cum.sums) |> head()
```

```
##  1  2  3  4  5  6 
##  1  1  1  1 72  1
```

```r
getQuantilepos_Vectorized_max_decoded <- Vectorize(vectorize.args = "cum.sum", function(q, cum.sum) {
  cum.sum <- decode(cum.sum) 
  c.max <- max(cum.sum) # Max is last element since x is a cumulative sums.
  treshold <- c.max * q
  which.max(cum.sum >= treshold)
})
getQuantilepos_Vectorized_max_decoded(.9, cum.sums) |> head()
```

```
##  1  2  3  4  5  6 
##  1  1  1  1 72  1
```

```r
getQuantilepos_Vectorized_max <- Vectorize(vectorize.args = "cum.sum", function(q, cum.sum) {
  c.max <- max(cum.sum) # Max is last element since x is a cumulative sums.
  treshold <- c.max * q
  which.max(cum.sum >= treshold)
})
getQuantilepos_Vectorized_max(.9, cum.sums) |> head()
```

```
##  1  2  3  4  5  6 
##  1  1  1  1 72  1
```

```r
getQuantilepos_lapply <- function(q, cum.sums, decode = FALSE) {
  if (isTRUE(decode)) cum.sums <- lapply(cum.sums, decode)
  l <- sapply(cum.sums, \(cum.sum) {
  c.max <- tail(cum.sum,1) # Max is last element since x is a cumulative sums.
  treshold <- c.max * q
  which.max(cum.sum >= treshold)})
  l
}
getQuantilepos_lapply(.9, cum.sums) |> head()
```

```
##  1  2  3  4  5  6 
##  1  1  1  1 72  1
```

```r
getQuantilepos_Rle_native <- function(q, cum.sums) {
  (cum.sums >= max(cum.sums) * q )|> sapply(which.max)
}
getQuantilepos_Rle_native(.9, cum.sums) |> head()
```

```
##  1  2  3  4  5  6 
##  1  1  1  1 72  1
```

```r
getQuantilepos_NL_native <- function(q, cum.sums) {
  nl <- as(cum.sums, "NumericList")
  min(which(nl >= max(nl) * q))
}
getQuantilepos_NL_native(.9, cum.sums) |> head()
```

```
##  1  2  3  4  5  6 
##  1  1  1  1 72  1
```

```r
getQuantilepos_ChatGPT  <- function(q, cum.sums, decode = FALSE) {
  if (isTRUE(decode)) cum.sums <- lapply(cum.sums, decode)
  find_quantile_and_index <- function(x, q) {
    # Interestingly ChatGPT's definition of quantiles appears more sound,
    # quantile_value <- quantile(x, q)
    # But the definition below is the one orgiginally used in CAGEr.
    quantile_value <- max(x) * q
    which(x >= quantile_value)[1]
  }
  sapply(cum.sums, find_quantile_and_index, q)
}
getQuantilepos_ChatGPT(.9, cum.sums) |> head()
```

```
##  1  2  3  4  5  6 
##  1  1  1  1 72  1
```

## Benchmark


```r
(microbench_out <- microbenchmark::microbenchmark(
  times = 100,
  Vectorized_tail_decoded = getQuantilepos_Vectorized_tail_decoded(.9, cum.sums),
  Vectorized_tail         = getQuantilepos_Vectorized_tail(.9, cum.sums),
  Vectorized_max_decoded  = getQuantilepos_Vectorized_max_decoded(.9, cum.sums),
  Vectorized_max          = getQuantilepos_Vectorized_max(.9, cum.sums),
  lapply                  = getQuantilepos_lapply(.9, cum.sums),
  lapply_decoded          = getQuantilepos_lapply(.9, cum.sums, decode = TRUE),
  Rle_native              = getQuantilepos_Rle_native(.9, cum.sums),
  NL_native               = getQuantilepos_NL_native(.9, cum.sums),
  ChatGPT                 = getQuantilepos_ChatGPT(.9, cum.sums),
  ChatGPT_decoded         = getQuantilepos_ChatGPT(.9, cum.sums, decode = TRUE)
))
```

```
## Unit: milliseconds
##                     expr         min          lq       mean     median
##  Vectorized_tail_decoded  226.829726  234.819918  286.97939  253.23864
##          Vectorized_tail  684.732639  735.672604  904.44242  827.81918
##   Vectorized_max_decoded  220.561623  226.360906  286.71244  250.03909
##           Vectorized_max  425.882113  451.098007  585.66493  511.43531
##                   lapply 1218.524937 1358.116061 1578.09167 1482.08125
##           lapply_decoded  679.861133  735.194260  879.32682  817.87322
##               Rle_native  675.975803  718.836427  838.95589  779.28237
##                NL_native    9.486656    9.942683   11.53794   10.42981
##                  ChatGPT  952.030342 1015.902702 1277.46268 1122.26911
##          ChatGPT_decoded  677.942011  742.759078  875.11836  816.32189
##          uq        max neval
##   294.79168  700.34834   100
##   896.60301 1709.49658   100
##   301.12672  892.87669   100
##   655.94249 1334.79107   100
##  1602.81593 3465.68446   100
##   912.35238 2064.05496   100
##   872.24652 2335.11442   100
##    11.83217   26.55179   100
##  1345.28727 7607.42296   100
##   890.72910 1632.77436   100
```

```r
# https://statisticsglobe.com/microbenchmark-package-r
ggplot(microbench_out, aes(x = time / 1e9, y = expr, color = expr)) +  # Plot performance comparison
  geom_boxplot() + 
  scale_x_log10("time (seconds)")
```

![](quantile_position_files/figure-html/benchmark-1.png)<!-- -->

## Result

The winner is 20 times faster!

```
getQuantilepos_NL_native <- function(q, cum.sums) {
  nl <- as(cum.sums, "NumericList")
  min(which(nl >= max(nl) * q))
}
```


```r
identical(
  getQuantilepos_Vectorized_tail_decoded(.9, cum.sums),
  getQuantilepos_NL_native(.9, cum.sums)
)
```

```
## [1] TRUE
```

## Session information


```r
sessionInfo()
```

```
## R version 4.3.0 (2023-04-21)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Debian GNU/Linux 12 (bookworm)
## 
## Matrix products: default
## BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.11.0 
## LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.11.0
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
## 
## time zone: Etc/UTC
## tzcode source: system (glibc)
## 
## attached base packages:
## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
## [8] base     
## 
## other attached packages:
##  [1] ggplot2_3.4.2               CAGEr_2.6.0                
##  [3] MultiAssayExperiment_1.26.0 SummarizedExperiment_1.30.2
##  [5] Biobase_2.60.0              GenomicRanges_1.52.0       
##  [7] GenomeInfoDb_1.36.1         IRanges_2.34.1             
##  [9] S4Vectors_0.38.1            BiocGenerics_0.46.0        
## [11] MatrixGenerics_1.12.2       matrixStats_1.0.0          
## 
## loaded via a namespace (and not attached):
##   [1] RColorBrewer_1.1-3        rstudioapi_0.14          
##   [3] jsonlite_1.8.5            magrittr_2.0.3           
##   [5] GenomicFeatures_1.52.1    farver_2.1.1             
##   [7] rmarkdown_2.22            BiocIO_1.10.0            
##   [9] zlibbioc_1.46.0           vctrs_0.6.2              
##  [11] memoise_2.0.1             Rsamtools_2.16.0         
##  [13] DelayedMatrixStats_1.22.1 RCurl_1.98-1.12          
##  [15] base64enc_0.1-3           htmltools_0.5.5          
##  [17] S4Arrays_1.0.4            progress_1.2.2           
##  [19] curl_5.0.1                Formula_1.2-5            
##  [21] sass_0.4.6                KernSmooth_2.23-20       
##  [23] bslib_0.5.0               htmlwidgets_1.6.2        
##  [25] plyr_1.8.8                Gviz_1.44.0              
##  [27] cachem_1.0.8              GenomicAlignments_1.36.0 
##  [29] lifecycle_1.0.3           pkgconfig_2.0.3          
##  [31] Matrix_1.5-3              R6_2.5.1                 
##  [33] fastmap_1.1.1             GenomeInfoDbData_1.2.10  
##  [35] digest_0.6.31             colorspace_2.1-0         
##  [37] AnnotationDbi_1.62.1      Hmisc_5.1-0              
##  [39] RSQLite_2.3.1             vegan_2.6-4              
##  [41] filelock_1.0.2            fansi_1.0.4              
##  [43] mgcv_1.8-41               httr_1.4.6               
##  [45] compiler_4.3.0            microbenchmark_1.4.10    
##  [47] withr_2.5.0               bit64_4.0.5              
##  [49] htmlTable_2.4.1           backports_1.4.1          
##  [51] CAGEfightR_1.20.0         BiocParallel_1.34.2      
##  [53] DBI_1.1.3                 highr_0.10               
##  [55] biomaRt_2.56.1            MASS_7.3-58.2            
##  [57] rappdirs_0.3.3            DelayedArray_0.26.5      
##  [59] rjson_0.2.21              permute_0.9-7            
##  [61] gtools_3.9.4              tools_4.3.0              
##  [63] foreign_0.8-84            nnet_7.3-18              
##  [65] glue_1.6.2                restfulr_0.0.15          
##  [67] nlme_3.1-162              stringdist_0.9.10        
##  [69] grid_4.3.0                checkmate_2.2.0          
##  [71] reshape2_1.4.4            cluster_2.1.4            
##  [73] generics_0.1.3            operator.tools_1.6.3     
##  [75] gtable_0.3.3              BSgenome_1.68.0          
##  [77] formula.tools_1.7.1       ensembldb_2.24.0         
##  [79] data.table_1.14.8         hms_1.1.3                
##  [81] xml2_1.3.4                utf8_1.2.3               
##  [83] XVector_0.40.0            pillar_1.9.0             
##  [85] stringr_1.5.0             splines_4.3.0            
##  [87] dplyr_1.1.2               BiocFileCache_2.8.0      
##  [89] lattice_0.20-45           rtracklayer_1.60.0       
##  [91] bit_4.0.5                 deldir_1.0-9             
##  [93] biovizBase_1.48.0         tidyselect_1.2.0         
##  [95] Biostrings_2.68.1         knitr_1.43               
##  [97] gridExtra_2.3             ProtGenerics_1.32.0      
##  [99] xfun_0.39                 stringi_1.7.12           
## [101] VGAM_1.1-8                lazyeval_0.2.2           
## [103] yaml_2.3.7                som_0.3-5.1              
## [105] evaluate_0.21             codetools_0.2-19         
## [107] interp_1.1-4              tibble_3.2.1             
## [109] cli_3.6.1                 rpart_4.1.19             
## [111] munsell_0.5.0             jquerylib_0.1.4          
## [113] dichromat_2.0-0.1         Rcpp_1.0.10              
## [115] dbplyr_2.3.2              png_0.1-8                
## [117] XML_3.99-0.14             parallel_4.3.0           
## [119] assertthat_0.2.1          blob_1.2.4               
## [121] prettyunits_1.1.1         latticeExtra_0.6-30      
## [123] jpeg_0.1-10               AnnotationFilter_1.24.0  
## [125] sparseMatrixStats_1.12.1  bitops_1.0-7             
## [127] VariantAnnotation_1.46.0  scales_1.2.1             
## [129] crayon_1.5.2              rlang_1.1.1              
## [131] KEGGREST_1.40.0
```

## ChatGPT

Charles  10:15 AM

> Good morning ChatGPT

APP  10:15 AM

> Good morning! How can I assist you today?

Charles   10:17 AM

> Using the R programming language, I want to find in a list of vectors of sorted numerical values the index where the value exceeds a given quantile.

APP  10:17 AM

> In R programming, you can use the `which` function in combination with the `quantile` function to find the index in a list of sorted numerical vectors where the value exceeds a given quantile. Here's an example:

```
r
# Assume 'v' is your sorted numerical vector
v <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)

# Define your quantile, for example, the 50th percentile (median)
q <- quantile(v, 0.5)

# Find the indices where the value exceeds the quantile
indices <- which(v > q)

# Print the indices
print(indices)
```

> This will give you the indices of the values in the vector 'v' that are greater than the 50th percentile. You can adjust the quantile value as needed.

Charles  10:18 AM

> Thanks, the vectors are in a list, and I need to run the function on each of them in an efficient way.  Can you write a vectorized function?

APP  10:18 AM

> In R programming, you can use the `lapply` function in combination with which and quantile to create a vectorized function that finds the index in each vector of a list `where` the value exceeds a given `quantile.` Here's an example:

```
r
# Define your list of vectors
list_of_vectors <- list(c(1,2,3,4,5), c(6,7,8,9,10), c(11,12,13,14,15))

# Define the quantile
q <- 0.75

# Define the function
find_index <- function(x, q) {
  return(which(x > quantile(x, q)))
}

# Apply the function to each vector in the list
result <- lapply(list_of_vectors, find_index, q)

# Print the result
print(result)
```

> In this example, `find_index` is a function that takes a vector and a quantile as input, and returns the indices of the vector where the value exceeds the given quantile. The `lapply` function applies this function to each vector in the list. The result is a list of vectors, where each vector contains the indices where the value exceeds the given quantile.

Charles  10:19 AM

> Thanks, we are almost there. I want the function to compute the quantile specifically on each vector.

APP  10:19 AM

> In R programming, you can use the `lapply()` function to apply a function to each element of a list. If you want to compute the quantile for each vector in the list, you can do it like this:

```
r
# Assuming 'my_list' is your list of vectors
result <- lapply(my_list, function(x) quantile(x, probs = c(0.25, 0.5, 0.75)))
```

> This will return a new list where each element is the 25th, 50th, and 75th percentile of the corresponding vector in the original list. You can modify the `probs` parameter to compute different quantiles.

Charles 10:21 AM

> Sorry, but what I mean is: I have a list of vectors, and for each vector of that list I need to find a given quantile, and the first index as which the value in the vector exceed the quantile.

APP  10:21 AM

> In R, you can use the `lapply` function to apply a function to each element in a list. Here's how you can find for each vector in a list, a given quantile, and the first index at which the value in the vector exceeds the quantile:

```
r
# Assume 'lst' is your list of vectors and 'q' is your quantile
lst <- list(c(1,2,3,4,5), c(6,7,8,9,10), c(11,12,13,14,15))
q <- 0.5

# Function to calculate quantile and find the first index exceeding the quantile
find_quantile_and_index <- function(x, q) {
  quantile_value <- quantile(x, q)
  index <- which(x > quantile_value)[1]
  return(list("quantile" = quantile_value, "index" = index))
}

# Apply the function to each vector in the list
result <- lapply(lst, find_quantile_and_index, q)

# Print the result
print(result)
```

> This will return a list of lists, where each sub-list contains the quantile and the first index at which the value in the vector exceeds the quantile for each vector in the original list.

Charles  10:29 AM

> Thanks, but it is too slow in my benchmark. Can you propose me a faster version ?
APP  10:29 AM

> To optimize your function in R, you can use the `lapply` function which is faster than loops. Here's an example:

```
r
# Assuming 'lst' is your list of vectors and 'q' is your quantile
results <- lapply(lst, function(x) {
  quantile_val <- quantile(x, q)
  index <- which(x > quantile_val)[1]
  return(list("quantile" = quantile_val, "index" = index))
})
```

> This function applies the anonymous function to each vector in the list. The function calculates the quantile, finds the first index where the value exceeds the quantile, and returns them as a list. The result is a list of these lists.

> Remember, R is vectorized which means it's more efficient to apply a function to a whole vector at once rather than element by element. Functions like `lapply` take advantage of this. Also, using `which` with the condition directly is faster than looping through the vector to find the index.
