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
```

The DataFrame encodes counts of CAGE tags mapped on the zebrafish genome.


```r
# DF <- CTSStagCountDF(exampleCAGEexp)
DF <- ZebrafishDevelopmentalCAGE::ZebrafishCAGE |> assay()
DF
```

```
## DataFrame with 3552111 rows and 12 columns
##         zf_unfertilized_egg zf_fertilized_egg zf_64cells zf_512cells zf_high
##                       <Rle>             <Rle>      <Rle>       <Rle>   <Rle>
## 1                         0                 2          0           0       0
## 2                         0                 0          1           0       0
## 3                         0                 0          0           0       0
## 4                         0                 0          0           0       0
## 5                         0                 0          0           0       0
## ...                     ...               ...        ...         ...     ...
## 3552107                   0                 0          0           1       0
## 3552108                   0                 0          0           0       0
## 3552109                   0                 0          0           0       0
## 3552110                   1                 0          0           0       0
## 3552111                   0                 0          0           0       0
##         zf_oblong zf_sphere_dome zf_30perc_dome zf_shield zf_14somites zf_prim6
##             <Rle>          <Rle>          <Rle>     <Rle>        <Rle>    <Rle>
## 1               0              0              0         0            0        0
## 2               0              0              0         0            0        0
## 3               0              0              0         1            0        0
## 4               0              0              0         0            0        0
## 5               0              0              0         0            0        1
## ...           ...            ...            ...       ...          ...      ...
## 3552107         0              0              0         0            0        0
## 3552108         0              0              0         0            1        0
## 3552109         0              1              0         0            0        0
## 3552110         0              0              0         0            0        0
## 3552111         0              0              0         0            1        0
##         zf_prim20
##             <Rle>
## 1               0
## 2               0
## 3               0
## 4               1
## 5               2
## ...           ...
## 3552107         0
## 3552108         0
## 3552109         0
## 3552110         0
## 3552111         0
```

```r
DF[[1]]
```

```
## integer-Rle of length 3552111 with 796039 runs
##   Lengths: 53  1 64  1  9  1  1  1  1  1  2 ...  3  1  1  3  1 33  1 78  1  1
##   Values :  0  1  0  1  0  1  0  1  0  5  0 ...  0  2  1  0  2  0  1  0  1  0
```

The functions:


```r
f_Reduce_native         <- function(DF) Reduce(`+`, DF, Rle(0L))
f_Reduce_native(DF)
```

```
## integer-Rle of length 3552111 with 2254482 runs
##   Lengths:  1  3  1  2  1  6  1  3  1  3  1 ...  3  1  2  1  3  1  2  1  1  7
##   Values :  2  1  3  1  5  1  6  1  2  1  3 ...  1  2  1  2  1  2  1  8  2  1
```

```r
f_decode_reduce         <- function(DF) Rle(Reduce(`+`, lapply(DF, decode), 0L))
f_decode_reduce(DF)
```

```
## integer-Rle of length 3552111 with 2254482 runs
##   Lengths:  1  3  1  2  1  6  1  3  1  3  1 ...  3  1  2  1  3  1  2  1  1  7
##   Values :  2  1  3  1  5  1  6  1  2  1  3 ...  1  2  1  2  1  2  1  8  2  1
```

```r
f_rowSums_as.data.frame <- function(DF) Rle(rowSums(as.data.frame(DF)))
f_rowSums_as.data.frame(DF)
```

```
## numeric-Rle of length 3552111 with 2254482 runs
##   Lengths:  1  3  1  2  1  6  1  3  1  3  1 ...  3  1  2  1  3  1  2  1  1  7
##   Values :  2  1  3  1  5  1  6  1  2  1  3 ...  1  2  1  2  1  2  1  8  2  1
```

```r
f_rowSums_decode        <- function(DF) Rle(rowSums(data.frame(lapply(DF, decode))))
f_rowSums_decode(DF)
```

```
## numeric-Rle of length 3552111 with 2254482 runs
##   Lengths:  1  3  1  2  1  6  1  3  1  3  1 ...  3  1  2  1  3  1  2  1  1  7
##   Values :  2  1  3  1  5  1  6  1  2  1  3 ...  1  2  1  2  1  2  1  8  2  1
```

```r
f_rowSums_delayedarray  <- function(DF) Rle(rowSums(DelayedArray::DelayedArray(DF)))
f_rowSums_delayedarray(DF)
```

```
## numeric-Rle of length 3552111 with 2254482 runs
##   Lengths:  1  3  1  2  1  6  1  3  1  3  1 ...  3  1  2  1  3  1  2  1  1  7
##   Values :  2  1  3  1  5  1  6  1  2  1  3 ...  1  2  1  2  1  2  1  8  2  1
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
##          f_Reduce_native(DF) 3046.0176 3066.7777 3211.4202 3115.1503 3230.7340
##          f_decode_reduce(DF)  404.7085  413.9336  504.3245  433.5773  476.5751
##  f_rowSums_as.data.frame(DF)  473.6294  493.5210  662.8725  545.0077  616.7744
##         f_rowSums_decode(DF)  463.1525  493.0847  668.2775  568.0069  637.9601
##   f_rowSums_delayedarray(DF) 2956.4122 3107.7231 3645.9180 3756.5327 3894.6892
##       max neval
##  3970.646   100
##  1157.323   100
##  1338.759   100
##  1333.968   100
##  5088.981   100
```

```r
# https://statisticsglobe.com/microbenchmark-package-r
ggplot(microbench_out, aes(x = time / 1e9, y = expr, color = expr)) +  # Plot performance comparison
  geom_boxplot() + 
  scale_x_log10("time (seconds)")
```

![](rowSums_on_Rle_DF_files/figure-html/benchmark-1.png)<!-- -->

## Result

The winner is:

```
f_decode_reduce <- function(DF) Rle(Reduce(`+`, lapply(DF, decode), 0L))
```

Explained with more verbose code:

```
f_decode_reduce <- function(DF) {
  list_of_numerical_vectors <- lapply(DF, decode)
  parallel_sum              <- Reduce(`+`, list_of_numerical_vectors, 0L)
  result                    <- Rle(parallel_sum)
  return(result)
}
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
##   [1] RColorBrewer_1.1-3                 rstudioapi_0.14                   
##   [3] jsonlite_1.8.5                     magrittr_2.0.3                    
##   [5] GenomicFeatures_1.52.1             farver_2.1.1                      
##   [7] rmarkdown_2.22                     BiocIO_1.10.0                     
##   [9] zlibbioc_1.46.0                    vctrs_0.6.2                       
##  [11] memoise_2.0.1                      Rsamtools_2.16.0                  
##  [13] DelayedMatrixStats_1.22.1          RCurl_1.98-1.12                   
##  [15] base64enc_0.1-3                    htmltools_0.5.5                   
##  [17] S4Arrays_1.0.4                     progress_1.2.2                    
##  [19] curl_5.0.1                         Formula_1.2-5                     
##  [21] sass_0.4.6                         KernSmooth_2.23-20                
##  [23] bslib_0.5.0                        htmlwidgets_1.6.2                 
##  [25] plyr_1.8.8                         Gviz_1.44.0                       
##  [27] cachem_1.0.8                       GenomicAlignments_1.36.0          
##  [29] lifecycle_1.0.3                    pkgconfig_2.0.3                   
##  [31] Matrix_1.5-3                       R6_2.5.1                          
##  [33] fastmap_1.1.1                      GenomeInfoDbData_1.2.10           
##  [35] digest_0.6.31                      colorspace_2.1-0                  
##  [37] AnnotationDbi_1.62.1               Hmisc_5.1-0                       
##  [39] RSQLite_2.3.1                      vegan_2.6-4                       
##  [41] filelock_1.0.2                     fansi_1.0.4                       
##  [43] mgcv_1.8-41                        httr_1.4.6                        
##  [45] compiler_4.3.0                     microbenchmark_1.4.10             
##  [47] withr_2.5.0                        bit64_4.0.5                       
##  [49] htmlTable_2.4.1                    backports_1.4.1                   
##  [51] CAGEfightR_1.20.0                  BiocParallel_1.34.2               
##  [53] DBI_1.1.3                          highr_0.10                        
##  [55] biomaRt_2.56.1                     MASS_7.3-58.2                     
##  [57] rappdirs_0.3.3                     DelayedArray_0.26.4               
##  [59] rjson_0.2.21                       permute_0.9-7                     
##  [61] gtools_3.9.4                       tools_4.3.0                       
##  [63] foreign_0.8-84                     nnet_7.3-18                       
##  [65] glue_1.6.2                         restfulr_0.0.15                   
##  [67] nlme_3.1-162                       stringdist_0.9.10                 
##  [69] grid_4.3.0                         checkmate_2.2.0                   
##  [71] reshape2_1.4.4                     cluster_2.1.4                     
##  [73] generics_0.1.3                     operator.tools_1.6.3              
##  [75] gtable_0.3.3                       BSgenome_1.68.0                   
##  [77] formula.tools_1.7.1                ensembldb_2.24.0                  
##  [79] ZebrafishDevelopmentalCAGE_0.999.0 data.table_1.14.8                 
##  [81] hms_1.1.3                          xml2_1.3.4                        
##  [83] utf8_1.2.3                         XVector_0.40.0                    
##  [85] pillar_1.9.0                       stringr_1.5.0                     
##  [87] splines_4.3.0                      dplyr_1.1.2                       
##  [89] BiocFileCache_2.8.0                lattice_0.20-45                   
##  [91] rtracklayer_1.60.0                 bit_4.0.5                         
##  [93] deldir_1.0-9                       biovizBase_1.48.0                 
##  [95] tidyselect_1.2.0                   Biostrings_2.68.1                 
##  [97] knitr_1.43                         gridExtra_2.3                     
##  [99] ProtGenerics_1.32.0                xfun_0.39                         
## [101] stringi_1.7.12                     VGAM_1.1-8                        
## [103] lazyeval_0.2.2                     yaml_2.3.7                        
## [105] som_0.3-5.1                        evaluate_0.21                     
## [107] codetools_0.2-19                   interp_1.1-4                      
## [109] tibble_3.2.1                       cli_3.6.1                         
## [111] rpart_4.1.19                       munsell_0.5.0                     
## [113] jquerylib_0.1.4                    dichromat_2.0-0.1                 
## [115] Rcpp_1.0.10                        dbplyr_2.3.2                      
## [117] png_0.1-8                          XML_3.99-0.14                     
## [119] parallel_4.3.0                     assertthat_0.2.1                  
## [121] blob_1.2.4                         prettyunits_1.1.1                 
## [123] latticeExtra_0.6-30                jpeg_0.1-10                       
## [125] AnnotationFilter_1.24.0            sparseMatrixStats_1.12.1          
## [127] bitops_1.0-7                       VariantAnnotation_1.46.0          
## [129] scales_1.2.1                       crayon_1.5.2                      
## [131] rlang_1.1.1                        KEGGREST_1.40.0
```

## Appendix: asking ChatGPT

After running this benchmark, I asked ChatGPT via Sidekick.  Here is our conversation:

Charles Plessy  10:05 PM

> Good evening, ChatGPT

Sidekick APP  10:05 PM

> Good evening! How can I assist you today?

Charles Plessy  10:06 PM

> Do you know the DataFrame class in the Bioconductor S4Vectors R package?

Sidekick APP  10:06 PM

> Yes, I'm familiar with it. The DataFrame class in the Bioconductor S4Vectors R package is an S4 class for representing a list of vectors of the same length. It's similar to a data frame in base R, but with more features. It's designed to handle genomic data, but can be used for any data that can be organized as a collection of equal-length vectors.

Charles Plessy  10:07 PM

> I use it to store numeric values compressed with the Rle class of the S4Vectors R package.

Sidekick APP  10:07 PM

> Sure, I'd be happy to help!

> The DataFrame class is a part of the S4Vectors package in R, which is a part of the Bioconductor project. It's designed to handle large amounts of data efficiently. It's similar to data.frames in base R, but with some additional features.

> The Rle (Run Length Encoding) class is also a part of the S4Vectors package. It's used to represent a sequence of repeated values in a more memory-efficient way. For example, instead of storing the sequence c(1, 1, 1, 2, 2, 3, 3, 3, 3), you could store it as an Rle object as (3 times 1, 2 times 2, 4 times 3).

> To store numeric values compressed with the Rle class in a DataFrame, you would first create the Rle object, and then include it in the DataFrame. 

> Here's a basic example:

```
R
# Load the S4Vectors package
library(S4Vectors)

# Create an Rle object
rle_object <- Rle(c(1, 1, 1, 2, 2, 3, 3, 3, 3))

# Create a DataFrame with the Rle object
df <- DataFrame(rle_column = rle_object)
```

> In this example, 'df' is a DataFrame that contains a single column ('rle_column') with the compressed numeric values.

> Remember that Rle objects can be more memory-efficient when you have long sequences of repeated values. If your data doesn't have many repeated values, you might not see much benefit from using Rle.

Charles Plessy  10:09 PM

> Can the base function rowSums work on a Bioconductor DataFrame of Rle-encoded values?

Sidekick APP  10:09 PM

> _Thinking... Thanks for waiting!_

Charles Plessy  10:20 PM

> I am waiting but please allow me to repeat my question: can the base function rowSums work on a Bioconductor DataFrame of Rle-encoded values?

Sidekick APP  10:20 PM

> No, the base function `rowSums` in R cannot directly work on a Bioconductor DataFrame of Rle-encoded values. The `rowSums` function is designed to work with numeric matrices or data frames. Rle-encoded values are a type of list, and `rowSums` cannot handle lists directly. You would need to first convert or extract the Rle-encoded values into a numeric matrix or data frame before using `rowSums`.

Charles Plessy  10:22 PM

> Please write for me an efficient function that computes row sums on a Bioconductor DataFrame of Rle-encoded values.

Sidekick APP  10:22 PM

> Sure, here's a simple function in R that computes the row sums on a DataFrame of Rle-encoded values. This function assumes that your DataFrame is named `df` and that it only contains Rle-encoded columns.

```
R
library(DelayedArray)

rowSumsRle <- function(df) {
  # Convert the DataFrame to a matrix
  mat <- as.matrix(df)
  
  # Apply the rowSums function
  result <- rowSums(mat)
  
  return(result)
}
You can use this function like this:
R
sums <- rowSumsRle(df)
```

> This will return a vector where each element is the sum of the corresponding row in the DataFrame. Please note that this function will only work correctly if all columns in the DataFrame are Rle-encoded. If there are any non-Rle-encoded columns, you will need to handle them separately.

Charles Plessy  10:26 PM

> Thank you, can you write one function that is more computationally efficient?

Sidekick APP  10:26 PM

> Sure, here's a simple example of how you might do this in R. This function uses the `rowSums` function, which is a highly optimized function for computing row sums in R. It first decodes the Rle-encoded values using the `runValue` function from the Bioconductor package, then computes the row sums.

```
R
library(Bioconductor)

computeRowSums <- function(df) {
  # Decode Rle-encoded values
  decoded_values <- lapply(df, function(x) { runValue(x) })
  
  # Convert list to data frame
  decoded_df <- do.call(cbind.data.frame, decoded_values)
  
  # Compute row sums
  row_sums <- rowSums(decoded_df)
  
  return(row_sums)
}
```

> Please note that this function assumes that all columns in the input data frame are Rle-encoded. If this is not the case, you may need to add some additional error checking to handle columns that are not Rle-encoded.

Charles Plessy  10:29 PM

> Thank you again. I would like to make a benchmark with different approaches. Can you propose a third alternative that uses a different algorithm^

Sidekick APP  10:29 PM

> Sure, here's an alternative function using the `rowSums` function in R. This function will first decode the Rle-encoded values, then compute the row sums. Please note that this function assumes that your DataFrame is named `df` and that all columns contain Rle-encoded values.

```
R
library(DelayedArray)

computeRowSums <- function(df) {
  # Decode Rle-encoded values
  decoded_df <- as.data.frame(lapply(df, function(x) unlist(runValue(x))))
  
  # Compute row sums
  row_sums <- rowSums(decoded_df)
  
  return(row_sums)
}

# Usage
row_sums <- computeRowSums(df)
```

> This function uses the `lapply` function to apply the decoding function to each column of the DataFrame. The `unlist` and `runValue` functions are used to decode the Rle-encoded values. The `rowSums` function is then used to compute the row sums. Please replace `df` with your actual DataFrame name.

### But the functions do not work...


```r
rowSumsRle <- function(df) {
  # Convert the DataFrame to a matrix
  mat <- as.matrix(df)
  
  # Apply the rowSums function
  result <- rowSums(mat)
  
  return(result)
}
```

```
> rowSumsRle(DF)

 Error in base::rowSums(x, na.rm = na.rm, dims = dims, ...) : 
  'x' must be an array of at least two dimensions
```


```r
computeRowSums <- function(df) {
  # Decode Rle-encoded values
  decoded_values <- lapply(df, function(x) { runValue(x) })
  
  # Convert list to data frame
  decoded_df <- do.call(cbind.data.frame, decoded_values)
  
  # Compute row sums
  row_sums <- rowSums(decoded_df)
  
  return(row_sums)
}
```

```
> computeRowSums(DF)
Error in data.frame(..., check.names = FALSE) : 
  arguments imply differing number of rows: 796039, 952755, 971432, 1077882, 890977, 930738, 933017, 926733, 579846, 1449206, 1308034, 1260173
```


```r
computeRowSums <- function(df) {
  # Decode Rle-encoded values
  decoded_df <- as.data.frame(lapply(df, function(x) unlist(runValue(x))))
  
  # Compute row sums
  row_sums <- rowSums(decoded_df)
  
  return(row_sums)
}
```

```
> computeRowSums(DF)

 Error in (function (..., row.names = NULL, check.rows = FALSE, check.names = TRUE,  : 
  arguments imply differing number of rows: 796039, 952755, 971432, 1077882, 890977, 930738, 933017, 926733, 579846, 1449206, 1308034, 1260173
```
