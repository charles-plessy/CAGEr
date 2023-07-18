# Changes in version 2.7.1

BUG FIXES

-   Correct quantile positions, which were shifted by one base.  This bug may
    have been introduced in version 1.22 or later.
-   Ensure cluster objects are properly sorted.  Fixes #79, introduced in
    version 2.6.0 and causing crashes or incorrect quantile calculations.
-   Apply fix for #77 (aggregateTagClusters losing TCs), which slipped out of
    2.6.0 because of Git branch mixup.
-   Fix _consensus cluster_ coordinates, where the `maxDist` padding was
    erroneously remaining in some parts of the computation.
-   Corrected on-the-fly cumulative sum computation for consensus clusters when
    `sample = NULL`.  The bug was causing incorrectly short quantile ranges.
-   Force the cluster names to stay sorted, to avoid a bug desynchronising
    quantile information and genome coordinates.
    
NEW FEATURES

-   Allow URLs to files in `getCTSS()` (Fixes #50).
-   Accelerated the computation of quantile position by ~20 times.
-   New `resetCAGEexp()` function.
-   New `flagByUpstreamSequences()` function.

# Changes in version 2.6.0

BUG FIXES

-   Re-enable the promoter shift functions (Fixes #4).
-   Compute dominant TSS information for *consensus clusters* (Fixes
    #1).

NEW FEATURES

-   Inter-quantile width is also computed when no sample is selected.
-   Enhancer detection using a wrapper to the _CAGEfighter_ package's
    function `quickEnhancers()`.

# No changes in version 2.4.0 (2/11/2022)

# Changes in version 2.2.0 (27/4/2022)

BUG FIXES

-   Restore the `correctSystematicG` option in `getCTSS()`. See #61.
-   Restore object class consistency in if / else statement in private
    function `bam2CTSS`. Fixes #49.
-   Restore proper CTSS conversion from BAM files (bug introduced in
    v1.34.0) while fixing issue #36.
-   Ensure Tag Clusters have a Seqinfo. Fixes #63.

# Changes in version 2.0.0 (27/10/2021)

BACKWARDS-INCOMPATIBLE CHANGES

-   The `CAGEset` class is removed.
-   Accessors using plain `data.frame` formats are removed.
-   Modifier functions return the object instead of silently modifying
    it in the global environment. Thus, you can use R pipes (`|>`).
-   Removed export function that unconditionally wrote files to the
    working directory, such as `exportCTSStoBedGraph`. As a replacement
    a new converter is provided, `exportToBrowserTrack`, that produces
    `UCSCData` objects that the user can then wite to files using the
    rtracklayer\_ package.
-   Removed the `extractExpressionClass` function, as the information is
    easily accessible from within the `CTSS` or `ConsensusClusters`
    objects.

NEW FEATURES

-   `CTSSnormalizedTpmGR` and `CTSStagCountGR` now accept `all` as a
    special sample name to return them all in a `GRangesList` object.
-   Plot interquantile width and expression clusters with *ggplot2*.

BUG FIXES

-   Corrected the `getExpressionProfiles` function to accept *CAGEexp*
    objects.
-   Updated the `exampleCAGEexp` object to contain expression classes.
-   Restore paraclu support for *CAGEexp* objects.
-   Corrected a bug in `aggregateTagClusters` that was causing
    mislabelling of tag clusters (PR#42).
-   Prevent `plotReverseCumulatives` from crashing when values are not
    in range. (PR#43).

# Changes in version 1.34.0 (20/5/2021)

BUG FIXES

-   Reform the *CTSS* class. New accessor: `CTSS()` (with no dot).
-   Correct a class error when loading BAM files. (Closes #36).
-   Use the *BSgenome* object from the main environment if available.

# Changes in version 1.32.0 (28/10/2020)

BUG FIXES

-   Update end coordinates before start coordinates in the function
    `.aggregateTagClustersGR()`. This should stop triggering
    *\"'width(x)' cannot contain negative integers\"* errors.
-   Correct `.make.consensus.clusters` internal function. This should
    stop triggering *\"Consensus clusters must not overlap with each
    other\"* errors.
-   Allow empty `CTSS.chr` objects.
-   Correct `plotInterquantileWidth()` to really use consensus clusters
    when passed the argument `clusters = "consensusClusters"`.
-   Fix failures on *CAGEexp* objects containing only one sample.

# Changes in version 1.30.0 (28/5/2020)

BUG FIXES

-   Correct usage of `&&` for vector comparison in pairs.DataFrame().
-   Adjust to latest *SummarizedExperiment* *MultiAssayExperiment*
    versions.

# Changes in version 1.28.0 (30/10/2019)

BACKWARDS-INCOMPATIBLE CHANGES

-   The *CTSS* objects are now based on the *UnstitchedGPos* class
    instead of the parent *GPos* class. Existing *CAGEexp* objects might
    give errors after package upgrade.

NEW FEATURES

-   Depends on *MultiAssayExperiment* instead of previously importing
    it. No need to load it by hand (nor *SummaryExperiment*) anymore.

BUG FIXES

-   Correct strand in `.remove.added.G` private function (PR#26).

# Changes in version 1.26.0 (3/5/2019)

DOCUMENTATION

-   Split the \"*CAGE resources*\" documentation into a separate
    vignette.

# Changes in version 1.24.0 (31/10/2018)

NEW FEATURES

-   `plotAnnot()`: the `group` argument now accepts formulas such as
    `~ a + b` to indicate names of metadata column that will be pasted
    together to form a new group factor.

BUG FIXES

-   Prevent `mergeSamples()` from producing `colData` that cause other
    functions to crash later when coercing to *data.frame*.
-   Repaired paraclu support for *CAGEset* objects.
-   `normalizeTagCount()` works again on *CAGEset* objects.
-   `consensusClustersGR()` reports expression score of the selected
    sample (instead of silently ignoring the `sample` argument and
    reporting expression sum on all the samples).

# Changes in version 1.22.0 (1/5/2018)

BACKWARDS-INCOMPATIBLE CHANGES

-   The plotting functions send their output to the graphical device
    instead of writing it to a file. This makes their use more
    consistent with most plotting functions in R.

NEW FEATURES

-   New *CAGEexp* class extending the *MultiAssayExperiment* class. It
    stores expression data more efficiently than *CAGEset*, and uses
    core Bioconductor types natively. For backwards compatibility it
    also support many of the original generic functions for *CAGEset*
    objects. Sample names in *CAGEexp* objects must by syntactically
    valid. Some methods are not yet available for *CAGEexp* object (see
    the vignette). Conversely some new functions will not be available
    for *CAGEset* objects.
-   New *CTSS* and *TagClusters* classes wrapping *GPos* and *GRanges*
    objects, for more type safety.
-   New functions for quality controls such as `plotAnnot()` or
    `hanabiPlot()`. See the *CAGEexp* vignette for details.
-   New `removeStrandInvaders()` function to count and remove
    strand-invasion artefacts (see Tang et al., 2013,
    <doi:10.1093/nar/gks112>).
-   Data export as *DESeqDataSet* object for *DESeq2* with the new
    `consensusClustersDESeq2()` and `GeneExpDESeq2()` functions.
-   New `bedctss` format to load the FANTOM5 and FANTOM6 CAGE data.
-   Multicore parallelisation with *BiocParallel* instead of *parallel*.
-   New function `sampleList()` to help looping on samples with
    `lapply()`.
-   New `plotCorrelation2()` function, faster than `plotCorrelation()`
    because it is plain black and white.
-   Multicore loading of CTSS data in *CAGEexp* objects.

OTHER CHANGES

-   Example data `exampleCAGEexp`, `exampleCAGEexp` and
    `exampleZv9_annot` are now is lazy-loaded.
-   `NULL` can be passed as genome name, to circumvent the requirement
    for a *BSgenome* object when actually not needing one.
-   In *CAGEexp* objects, expression quantile positions are given
    relative to the cluster start site.
-   For performance reasons, the positions of a quantile Q is now
    calculated as the position of the first base where cumulative
    expression is higher or equal to Q% of the total expression of a
    cluster.
-   `plotReverseCumulatives()`:
    -   fitInRange accepts a value of `NULL` to turn off power law
        fitting.
    -   New `legend` argument to remove legend when set to `FALSE`.
    -   Axis range and labels can be modified with xlab/ylab and
        xlim/ylim.
    -   Ladders steps start on the values instead of being centered on.
-   `setColors()`: allow lowercase in color names.

BUG FIXES

-   Corrected a bug that was crashing *CAGEset* objects when loading
    more than one BAM file.
-   Load BAM as gapped alignment with `readGAlignments()` instead of
    `scanBam()`. Without this correction, TSS position on minus strand
    is incorrect in case of indels in the read.

DOCUMENTATION UPDATES

-   Roxygen is used to generate the manual pages.
-   A new HTML vignette describes the *CAGEexp* class.

# Changes in version 1.18.1

-   `getCTSS()` can now load data from large BAM files. Before this
    update, there were "negative length vectors" errors when filtering
    out low quality reads from datasets with multiple dozens of millions
    of sequences.

# Changes in version 1.18.0 (25/4/2017)

-   Added Charles Plessy as co-maintainer.
-   Remove warning by replacing deprecated `ignoreSelf()` with
    `drop.self()`.

# Changes in version 1.13.4 (27/04/2016)

-   26/04/2016 updated the vignette.
-   26/04/2016 added support for BED file CAGE tags input - included the
    new accepted value for `inputFilesType="bed"`.
-   24/04/2016 added support for paired-end BAM file CAGE data input
    (CAGEscan). Included the new accepted value for
    `inputFilesType="bamPairedEnd"`.

# Changes in version 1.13.3 (18/04/2016)

-   18/04/2016 updated the vignette.
-   14/04/2016 modified the `consensusClusters()` function to allow
    extracting sample-specific information, such as CAGE signal,
    position of the dominant TSS and interquantile width in a specific
    sample (i.e. `consensusClusters()` now behaves analogously to
    `tagClusters()` function).
-   10/03/2016 fixed the bandwidth of the smoothing kernel for plotting
    correlation density plots to always scale with the number of bins.
-   10/03/2016 added new argument to `plotCorrelation()` function that
    allows controlling figure size in pixels.
-   03/11/2015 fixed the format of exporting bed12 files for tag
    clusters, which contained overlapping blocks per record and caused
    problems when trying to convert bed to bigBed.
-   03/11/2015 added an option of exporting CTSS bedGraph files to
    BigWig format - added a new argument to `exportCTSStoBedGraph()`
    function.

# No changes in version 1.13.0 (15/10/2015)

# Changes in version 1.10.2 (11/10/2015)

-   11/10/2015 changed to selective import from *data.table* package to
    resolve the conflict regarding the shift function.

# Changes in version 1.10.1 (17/09/2015)

-   15/05/2015 updated the vignette.
-   15/05/2015 enabled changing `xlim` and `ylim` when plotting
    interquantile width distribution.
-   13/05/2015 added `mergeCAGEsets()` function for merging two
    *CAGEset* objects.
-   11/05/2015 defined `as` method for coercing a data.frame into a
    *CAGEset* object.

# No changes in version 1.11.0 (17/04/2015)

# Changes in version 1.9.3 (06/02/2015)

-   06/02/2015 *CAGEr* package has now been published in NAR! If you
    find this package useful in your work, please cite: Haberle et al.
    (2015) CAGEr: precise TSS data retrieval and high-resolution
    promoterome mining for integrative analyses. Nucleic Acids Research
    43(8):e51, <doi:10.1093/nar/gkv054>. Check `citation("CAGEr")` for
    more details.
-   13/01/2015 fixed the error message in `tagClusters()` function for
    the case no `qLow` and `qUp` values are specified.

# Changes in version 1.9.2 (03/12/2014)

-   03/12/2014 updated the vignette - fixed the references.

# Changes in version 1.9.1 (01/12/2014)

-   01/12/2014 fixed the bug in importing data from ENCODEprojectCAGE
    data package using `importPublicData()` function when importing from
    multiple datasets that contain groups with the same name.

# No changes in version 1.9.0 (14/10/2014)

# Changes in version 1.7.2 (01/06/2014)

-   01/06/2014 replaced the .RData files in the /data folder with the
    more compressed ones.

# Changes in version 1.7.1 (30/05/2014)

-   30/05/2014 updated the vignette with the description of available
    public CAGE data resources and how to use them in *CAGEr*.
-   25/05/2014 added a list of available human and mouse FANTOM5 samples
    as a dataset that can be loaded via call to
    .`data(FANTOM5humanSamples)` or `data(FANTOM5mouseSamples)`.
-   25/05/2014 modified the `importPublicData()` function so that it
    supports importing FANTOM5 data for human and mouse directly from
    FANTOM web resource, as well as data from several CAGE data
    packages.
-   23/05/2014 optimized the creation of CTSS table across samples by
    performing merge on data.tables rather than *data.frames* - reading
    the data into *CAGEset* object from BAM or CTSS files is now \~10×
    faster.

# No changes in version 1.7.0 (14/04/2014)

# Changes in version 1.5.5 (01/04/2014)

-   13/02/2014 fixed bug in function for plotting interquantile width
    (in cases of multiple pairs of low and high quantiles present it was
    sometime fetching multiple columns and stopping with an error).
-   12/02/2014 fix typo in Details section of man page for `exportToBed`
    function.
-   12/02/2014 fixed the bug in function for systematic correction of G
    addition bias (it was stopping with an error in a rare case when all
    reads beginning with a G and mapping to the genome are derived by
    trimming the G mismatch from upstream position).

# Changes in version 1.5.3 (10/11/2013)

-   08/11/2013 fixed the bug in `importPublicData` function that caused
    problems when importing only one sample.
-   08/11/2013 fixed error in BibTex file.

# Changes in version 1.5.2 (06/11/2013)

-   06/11/2013 updated vignette and NEWS file.

# Changes in version 1.5.1 (02/11/2013)

-   02/11/2013 moved vignette source to /vignettes folder.

# No changes in version 1.5.0 (15/10/2013)

# Changes in version 1.3.10 (24/09/2013)

-   24/09/2013 optimized calculation of K-S p-values for shifting
    promoters using vectorized functions.
-   23/09/2013 removed a .C call to the pKS2 function from stats package
    and replaced it by an R implementation of pKS2.

# Changes in version 1.2.9 (23/08/2013)

-   23/08/2013 added `consensusClustersTpm()` function for retrieving
    matrix with TPM values for consensus clusters across all samples.

# Changes in version 1.2.8 (23/08/2013)

-   22/08/2013 modified all functions to work on data from only one
    strand (some were previously giving error in such case).
-   21/08/2013 fixed the bug in calling CTSSs from reads on minus strand
    when not correcting for G addition, i.e. when r`emoveFirstG=FALSE`.

# Changes in version 1.2.7 (08/08/2013)

-   08/08/2013 added new argument to `aggregateTagClusters()` function
    (`excludeSignalBelowThreshold`) that controls which TCs will
    contribute to the total CAGE signal of the consensus cluster - only
    the ones above the threshold that are used initially to set the
    boundaries of the consensus cluster, or all TCs that overlap the
    resulting consensus clusters.
-   07/08/2013 enabled retrieving interquantile width for tag clusters -
    added new returnInterquantileWidth argument to `tagClusters()`
    function.
-   06/07/2013 enabled importing entire CTSS table with multiple columns
    of tag counts representing multiple samples/experiments - allowed
    setting imputFilesType to "CTSStable" when creating CAGEset object
    and extended `getCTSS` function accordingly to read the entire table
    at once and fill the object.

# Changes in version 1.2.6 (23/07/2013)

-   22/07/2013 replaced *multicore* package dependency with *parallel*
    package (suggested by Bioconductor core team).
-   21/07/2013 fixed the bug in using loaded genome when correcting for
    G nucleotide bias.
-   20/07/2013 fixed plotReverseCumulatives to work with only one sample
    in the *CAGEset.*
-   19/07/2013 changed naming of correlation plot files to distinguish
    CTSS from consensus clusters plot.

# Changes in version 1.2.5 (28/06/2013)

-   28/06/2013 updated the vignette.

# Changes in version 1.2.4 (16/06/2013)

-   16/06/2013 allowed plotting reverse cumulatives for normalized tag
    counts - added new argument to `plotReverseCumulatives()` function.
-   15/06/2013 implemented calculating suggested referent power-law
    distribution and visualizing the distribution and its parameters on
    reverse cumulatives plot when calling `plotReverseCumulatives()`
    function.
-   14/06/2013 allowed assigning user-specified colors to samples to be
    used in visualizations - added new `setColors()` function.
-   13/06/2013 allowed plotting scatter plots and correlation for both
    individual TSSs (raw or normalized values) and consensus clusters -
    added new argument to `plotCorrelation()` function.

# Changes in version 1.2.3 (12/06/2013)

-   12/06/2013 implemented plotting pairwise scatter plots of CAGE tag
    count and calculating correlation between samples - added new
    `plotCorrelation()` function.

# Changes in version 1.2.2 (30/05/2013)

-   30/05/2013 implemented the algorithm for correcting 'G' nucleotide
    addition bias to CAGE tags described in Carninci et al., Nature
    Genetics 2006 - added a new option `correctSystematicG` to getCTSS
    function.

# Changes in version 1.2.1 (28/05/2013)

-   28/05/2013 fixed the bug in reporting number of CTSSs in cluster
    when clustering CTSSs using custom clusters (for empty clusters it
    was returning that the number of CTSSs is 1, although the signal is
    0)  
-   28/05/2013 fixed the bug in clustering CTSSs using custom clusters
    when run with `multicore=FALSE` (it was returning to many clusters -
    correct specified region but all of the regions on each chromosome,
    not only the specified one).
-   22/05/2013 fixed the column classes in CTSS data.frame when reading
    the data from bam files (to match the classes when reading from ctss
    files).
-   21/05/2013 updated the vignette with explanation on how to use
    custom build genomes.
-   21/05/2013 allowed usage of custom build *BSgenome* packages by
    removing the check whether the specified genome is present in the
    `available.genomes()` from *BSgenome* package.

# Changes in version 1.2.0 (02/05/2013)

-   02/05/2013 updated the vignette.
-   02/05/2013 added a new feature in the `plotReverseCumulatives`
    function, so that the slope of the fitted power-law distribution in
    the user selected range of values is shown on the plot for each
    sample (helps to choose appropriate alpha parameter for
    normalization).
-   02/05/2013 fixed a bug in selecting the range of tag count values
    for fitting power-law distribution in normalization.
-   01/05/2013 added sample labels as names to library sizes vector.
-   29/04/2013 implemented statistical testing (Kolmogorov-Smirnov test)
    for differential promoter usage based on cumulative sums of CAGE
    signal along the promoters (implemented within the `scoreShift`
    function).
-   20/04/2013 added pass-through option `none` for `method` parameter
    in `normalizeTagCount` function to enable using raw tag counts in
    downstream steps (CTSS clustering, promoter width, etc.).
-   19/04/2013 optimized `scoreShift` function to extract and process
    only the cumulative sums for samples being compared.
-   19/04/2013 replaced lapply with a for loop in the `scoreShift`
    function to avoid invoking `multicore` within `lapply`.
-   18/04/2013 fixed wrong error message in `plotExpressionProfiles`
    function that notifies about accepted values for `what` parameter.
-   27/03/2013 added *data.table* to the list of dependencies and
    optimized various parts of the code to use *data.table* instead of
    *data.frame*.
