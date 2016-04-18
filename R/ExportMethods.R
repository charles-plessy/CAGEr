###########################################################
# Functions for exporting results (graphical and textual)
#


setGeneric(
name="plotReverseCumulatives",
def=function(object, values = "raw", fitInRange = c(10, 1000), onePlot = FALSE){
	standardGeneric("plotReverseCumulatives")
}
)

setMethod("plotReverseCumulatives",
signature(object = "CAGEset"),
function (object, values = "raw", fitInRange = c(10, 1000), onePlot = FALSE){
		
	sample.labels <- sampleLabels(object)
	if(values == "raw"){
		tag.count <- object@tagCountMatrix
	}else if(values == "normalized"){
		tag.count <- object@normalizedTpmMatrix
	}else{
		stop("'values' parameter must be one of the (\"raw\", \"normalized\")")
	}
	
	pdf(file = paste("CTSS_reverse_cumulatives_", values, "_all_samples.pdf", sep = ""), width = 8, height = 8, onefile = T, bg = "transparent", family = "Helvetica", fonts = NULL)
	par(mar = c(5,5,5,2))
	cols <- names(sample.labels)
	
	if(values == "raw"){
		fit.coefs.m <- apply(tag.count, 2, function(x) {.fit.power.law.to.reverse.cumulative(values = as.integer(x), val.range = fitInRange)})
		fit.slopes <- fit.coefs.m[1,]
		names(fit.slopes) <- sample.labels
		reference.slope <- min(median(fit.slopes), -1.05)
		library.sizes <- librarySizes(object)
		reference.library.size <- 10^floor(log10(median(library.sizes)))
#reference.intercept <- log(reference.library.size/zeta(-1*reference.slope))  # intercept on natural logarithm scale
		reference.intercept <- log10(reference.library.size/zeta(-1*reference.slope))  # intercept on log10 scale used for plotting with abline
	}else if(values == "normalized"){
#		fit.coefs.m <- apply(tag.count, 2, function(x) {.fit.power.law.to.reverse.cumulative(values = x, val.range = fitInRange)})
	}
	
	if(onePlot == TRUE){
		vals <- tag.count[, 1]		
		if(values == "raw"){
			.plotReverseCumulative(values = as.integer(vals), col = cols[1], title = "All samples")
			if(length(sample.labels) > 1){
				sapply(c(2:length(sample.labels)), function(x) {vals <- as.integer(tag.count[, sample.labels[x]]); .plotReverseCumulative(values = vals, col = cols[x], add = TRUE)})
			}
			abline(v = fitInRange, lty = "dotted")
			abline(a = reference.intercept, b = reference.slope, col = "#7F7F7F7F", lty = "longdash")
			legend("topright", legend = paste("(", formatC(-1*fit.slopes, format = "f", digits = 2), ") ", sample.labels, sep = ""), bty = "n", col = cols, text.col = cols, lwd = 2, cex = 1.3, y.intersp = 1.2)
			legend("bottomleft", legend = c("Ref. distribution:", paste(" alpha = ", sprintf("%.2f", -1*reference.slope), sep = ""), paste(" T = ", reference.library.size, sep = "")), bty = "n", col = NA, text.col = "#7F7F7F", cex = 1.3, y.intersp = 1.2)
		}else if(values == "normalized"){
			.plotReverseCumulative(values = vals, col = cols[1], title = "All samples")
			if(length(sample.labels) > 1){			
				sapply(c(2:length(sample.labels)), function(x) {vals <- tag.count[, sample.labels[x]]; .plotReverseCumulative(values = vals, col = cols[x], add = TRUE)})
			}
			legend("topright", legend = sample.labels, bty = "n", col = cols, text.col = cols, lwd = 2, cex = 1.3, y.intersp = 1.2)
		}
	}else{
		if(values == "raw"){
			sapply(sample.labels, function(x) {vals <- as.integer(tag.count[, x]); .plotReverseCumulative(values = vals, col = cols[which(sample.labels == x)], title = x, col.title = cols[which(sample.labels == x)]); abline(v = fitInRange, lty = "dotted"); abline(a = reference.intercept, b = reference.slope, col = "#7F7F7F7F", lty = "longdash"); text(min(fitInRange), 10^6, labels = paste(" alpha =", formatC(-1*fit.slopes[x], format = "f", digits = 2), sep = " "), adj = c(0,1), col = cols[which(sample.labels == x)], cex = 1.3); legend("bottomleft", legend = c("Ref. distribution:", paste(" alpha = ", sprintf("%.2f", -1*reference.slope), sep = ""), paste(" T = ", reference.library.size, sep = "")), bty = "n", col = NA, text.col = "#7F7F7F", cex = 1.3, y.intersp = 1.2)})
		}else if(values == "normalized"){
			sapply(sample.labels, function(x) {vals <- tag.count[, x]; .plotReverseCumulative(values = vals, col = cols[which(sample.labels == x)], title = x, col.title = cols[which(sample.labels == x)])})
		}
	}
	dev.off()
	message("\nFile 'CTSS_reverse_cumulatives_", values, "_all_samples.pdf' has been created in your working directory (", getwd(), ")")
	
}
)


setGeneric(
name="exportCTSStoBedGraph",
def=function(object, values = "normalized", format = "BigWig", oneFile = TRUE){
	standardGeneric("exportCTSStoBedGraph")
}
)

setMethod("exportCTSStoBedGraph",
signature(object = "CAGEset"),
function (object, values = "normalized", format = "BigWig", oneFile = TRUE){
		
	sample.labels <- sampleLabels(object)
    
    if(format == "BigWig"){
        reference.genome <- genomeName(object)
        if(reference.genome %in% rownames(installed.packages()) == FALSE){
            stop("Reference genome is not installed! Please install required BSgenome package before exporting to BigWig files!")
        }else if(!paste("package:", reference.genome, sep = "") %in% search()){
            stop("Reference genome is not loaded! Load the genome by calling 'library(", reference.genome, ")'")
        }else{
            genome <- get(ls(paste("package:", reference.genome, sep="")))
        }
    }
    
	if(values == "raw") {
		data <- cbind(CTSScoordinates(object), object@tagCountMatrix)
        if(format == "BigWig"){
            .export.bw.all(data = data, sample.labels = sample.labels, v = values, genome = genome)
        }else if(format == "bedGraph"){
            .export.bedgraph.all(data = data, sample.labels = sample.labels, v = values, oneFile = oneFile)
        }else{
            stop("'format' parameter must be one of the (\"BigWig\", \"bedGraph\")")
        }
	}else if(values == "normalized"){
        data <- cbind(CTSScoordinates(object), as.data.frame(object@normalizedTpmMatrix))
        if(format == "BigWig"){
            .export.bw.all(data = data, sample.labels = sample.labels, v = values, genome = genome)
        }else if(format == "bedGraph"){
            .export.bedgraph.all(data = data, sample.labels = sample.labels, v = values, oneFile = oneFile)
        }else{
            stop("'format' parameter must be one of the (\"BigWig\", \"bedGraph\")")
        }
	}else{
		stop("'values' parameter must be one of the (\"raw\", \"normalized\")")
	}
	
    message("\n", format, " file(s) for CTSS ", values, " counts have been created in your working directory (", getwd(), ")")
    
	invisible(1)
	
}
)


setGeneric(
name="plotInterquantileWidth",
def=function(object, clusters, tpmThreshold = 5, qLow = 0.1, qUp = 0.9, xlim = c(0,150), ...){
	standardGeneric("plotInterquantileWidth")
}
)

setMethod("plotInterquantileWidth",
signature(object = "CAGEset"),
function (object, clusters, tpmThreshold, qLow, qUp, xlim = c(0,150), ...){
	
	sample.labels <- sampleLabels(object)
	cols <- names(sample.labels)
	
	if(clusters == "tagClusters"){	
		
		if(length(object@tagClustersQuantileLow)>0 & length(object@tagClustersQuantileUp)>0) {
			if(!(paste("q_", qLow, sep = "") %in% colnames(object@tagClustersQuantileLow[[1]]) & paste("q_", qUp, sep = "") %in% colnames(object@tagClustersQuantileUp[[1]]))){
				stop("No data for given quantile positions! Run 'quantilePositions()' function for desired quantiles first!")
			}
		}else{
			stop("No data for given quantile positions! Run 'quantilePositions()' function for desired quantiles first!")
		}
		
		filename <- "TC"
		q.low <- object@tagClustersQuantileLow
		q.up <- object@tagClustersQuantileUp
		idx.list <- lapply(as.list(sample.labels), function(x) {
			   
								cl <- tagClusters(object, sample = x)
								idx <- cl$tpm >= tpmThreshold
								return(idx)
						   
							}
						   )
	
	}else if (clusters == "consensusClusters"){
		
		if(length(object@consensusClustersQuantileLow)>0 & length(object@consensusClustersQuantileUp)>0) {
			
			if(!(paste("q_", qLow, sep = "") %in% colnames(object@consensusClustersQuantileLow[[1]]) & paste("q_", qUp, sep = "") %in% colnames(object@consensusClustersQuantileUp[[1]]))){
				stop("No data for given quantile positions! Run 'quantilePositions()' function for desired quantiles first!")
			}
		}else{
			stop("No data for given quantile positions! Run 'quantilePositions()' function for desired quantiles first!")
		}
		
		filename <- "consensusClusters_interquantile_width_all_samples.pdf"
		q.low <- object@consensusClustersQuantileLow
		q.up <- object@consensusClustersQuantileUp
		cl <- object@consensusClustersTpmMatrix
		idx.list <- lapply(as.list(sample.labels), function(x) {idx <- cl[,x][cl[,x] > 0] >= tpmThreshold})		
		
	}else{
		stop("'clusters' parameter must be one of the (\"tagClusters\", \"consensusClusters\")")
	}

	names(idx.list) <- sample.labels
	
	pdf(file = paste(clusters, "_interquantile_width_all_samples.pdf", sep = ""), width = 8, height = 8, onefile = T, bg = "transparent", family = "Helvetica", fonts = NULL)
	par(mar = c(5,5,5,1))
	sapply(sample.labels, function(x) {
		   
		   q.low.s <- q.low[[x]]
		   q.low.s <- as.integer(q.low.s[, which(colnames(q.low.s) == paste("q_", qLow, sep = ""))])
		   q.up.s <- q.up[[x]]
		   q.up.s <- as.integer(q.up.s[, which(colnames(q.up.s) == paste("q_", qUp, sep = ""))])
		   width <- q.up.s[idx.list[[x]]] - q.low.s[idx.list[[x]]] + 1
		   h <- hist(width, breaks = round(max(width)/2), plot = F)
		   h$counts <- h$counts/sum(h$counts)
		   col <- as.integer(col2rgb(cols[which(sample.labels == x)]))/255
		   plot(h, xlim = xlim, main = x, xlab = paste(clusters, " interquantile width q", qLow, "-q", qUp, " (bp)", sep = ""), ylab = "relative frequency", col = rgb(col[1], col[2], col[3], 0.5), border = cols[which(sample.labels == x)], cex.axis = 1.8, cex.lab = 1.8, cex.main = 2.5, col.main = cols[which(sample.labels == x)], ...)
		   
		   }
		   )
	dev.off()
	message("\nFile '", clusters, "_interquantile_width_all_samples.pdf' has been created in your working directory (", getwd(), ")")
	
}
)


setGeneric(
name="plotExpressionProfiles",
def=function(object, what){
	standardGeneric("plotExpressionProfiles")
}
)

setMethod("plotExpressionProfiles",
signature(object = "CAGEset"),
function (object, what){
	
	if(what == "CTSS") {
		
		cl <- object@CTSSexpressionClasses
		if(length(cl)>0){
			tpm.mx <- object@normalizedTpmMatrix
			tpm.mx <- tpm.mx[as.integer(names(cl)),]
			cl.method <- object@CTSSexpressionClusteringMethod
		}else{
			stop("No CTSS expression profiling has been done yet! Run getExpressionProfiles function first!")
		}
		
	}else if(what == "consensusClusters"){
		
		cl <- object@consensusClustersExpressionClasses
		if(length(cl)>0){
			tpm.mx <- object@consensusClustersTpmMatrix
			tpm.mx <- tpm.mx[as.integer(names(cl)),]
			cl.method <- object@consensusClustersExpressionClusteringMethod
		}else{
			stop("No consensusClusters expression profiling has been done yet! Run getExpressionProfiles function first!")
		}
		
	}else{
		stop("'what' parameter must be one of the (\"CTSS\", \"consensusClusters\")")
	}
	
	l <- .extract.cluster.info(tpm.mx, cl)
	cl <- l[[1]]
	m <- l[[2]]
	file_name = paste(what, "_expression_profiles.pdf", sep = "")
	pdf(file = file_name, height = 5.5 * (max(cl[,2]) + 1) + 3, width = 6 * (max(cl[,1]) + 1))
	par(omi = c(3,0.5,0.5,0.5), mfrow = c(max(cl[,2]) + 1, max(cl[,1]) + 1))
	suppressWarnings(.plot.clusters.beanplots(value.matrix = m, cl = cl, cl.method = cl.method, dim.som.x = max(cl[,1]) + 1, dim.som.y = max(cl[,2]) + 1, cex.axis = 3, las = 2))
	dev.off()
	message("\nFile '", file_name, "' has been created in your working directory (", getwd(), ")")	
	
}
)


setGeneric(
name="exportToBed",
def=function(object, what, qLow = NULL, qUp = NULL, colorByExpressionProfile = FALSE, oneFile = TRUE){
	standardGeneric("exportToBed")
}
)


setMethod("exportToBed",
signature(object = "CAGEset"),
function (object, what, qLow = NULL, qUp = NULL, colorByExpressionProfile = FALSE, oneFile = TRUE){
	
	sample.labels <- sampleLabels(object)

	if(what == "CTSS") {
		
		oneFile <- TRUE
		use.blocks <- F
		ctss <- object@CTSScoordinates
		#filtered_ctss <- object@filteredCTSSidx

		if(colorByExpressionProfile == TRUE){
			cl <- object@CTSSexpressionClasses
			n <- names(cl)
			cl <- .extract.cluster.info(cl = cl)
			cl <- data.frame(ctss = n, x.cor = cl[,1], y.cor = cl[,2])		
			ctss <- merge(cl, ctss, by.x = "ctss", by.y = 0, all.x = T, all.y = F)
			ctss <- data.frame(ctss = ctss$ctss, chr = ctss$chr, start = ctss$pos-1, end = ctss$pos, strand = ctss$strand, x.cor = ctss$x.cor, y.cor = ctss$y.cor)
			track.file <- "CTSS.colored.by.expression.profile.bed"
			track.names <- list("CTSS (colored by expression profile)")
			clustering.method <- object@CTSSexpressionClusteringMethod
		}else{
			ctss <- data.frame(chr = ctss$chr, start = ctss$pos-1, end = ctss$pos, strand = ctss$strand)
			track.file <- "CTSS.pooled.samples.bed"
			track.names <- list("CTSS (pooled samples)")
			#filtered_cols <- c("TRUE" = c("0,0,0"), "FALSE" = c("127,127,127"))
			#cols = list(filtered_cols[as.character(filtered_ctss)])
			cols = list(c("0,0,0"))
		}
		clusters.q.list = list(ctss)
		
	}else if(what == "tagClusters") {
		
		colorByExpressionProfile <- FALSE
		
		if(length(qLow) > 0 & (length(object@tagClustersQuantileLow)>0 & length(object@tagClustersQuantileLow)>0)) {
		
		if(paste("q_", qLow, sep = "") %in% colnames(object@tagClustersQuantileLow[[1]]) & paste("q_", qUp, sep = "") %in% colnames(object@tagClustersQuantileUp[[1]])){
		
		use.blocks <- T
		q.low <- object@tagClustersQuantileLow
		q.low <- lapply(q.low, function(x) {colnames(x)[2:ncol(x)] <- paste("qLow_", do.call(rbind, strsplit(colnames(x)[2:ncol(x)], split = "_", fixed = T))[,2], sep = ""); return(x)})
		q.up <- object@tagClustersQuantileUp
		q.up <- lapply(q.up, function(x) {colnames(x)[2:ncol(x)] <- paste("qUp_", do.call(rbind, strsplit(colnames(x)[2:ncol(x)], split = "_", fixed = T))[,2], sep = ""); return(x)})
		q <- lapply(as.list(1:length(sample.labels)), function(x) {merge(q.low[[x]], q.up[[x]], by.x = "cluster", by.y = "cluster")})
		names(q) <- sample.labels
		clusters <- lapply(as.list(sample.labels), function(x) {tagClusters(object, sample = x)})
		names(clusters) <- sample.labels
		clusters.q.list <- lapply(as.list(sample.labels), function(x) {merge(clusters[[x]], q[[x]], by.x = "cluster", by.y = "cluster")})
		track.names <- paste(sample.labels, paste(" (tag clusters (TC) q(", qLow, ")-q(",qUp,"))", sep = ""), sep = "")
		r <- paste(".qLow", qLow, "_qUp", qUp, sep = "")
			
		}else{
			stop("No data for given quantile positions! Run 'quantilePositions()' function for desired quantiles first, or omit 'qLow' and 'qUp' parameters to use start and end coordinates instead!")
		}
		}else if(length(qLow) > 0){
			stop("No data for given quantile positions! Run 'quantilePositions()' function for desired quantiles first, or omit 'qLow' and 'qUp' parameters to use start and end coordinates instead!")
		}else{
			use.blocks <- F
			clusters.q.list <- lapply(as.list(sample.labels), function(x) {tagClusters(object, sample = x)})			
			track.names <- paste(sample.labels, paste(" (tag clusters (TC))", sep = ""), sep = "")
			r <- ""
		}
		
		names(clusters.q.list) <- sample.labels
		itemRgb = FALSE
		cols <- names(sample.labels)
		cols <- as.list(apply(sapply(cols, function(x) {as.integer(col2rgb(x))}), 2, function(y) {paste(y, collapse = ",")}))
		names(cols) <- sample.labels

		if(oneFile){
			track.file <- rep(paste("All.samples.tagClusters", r, ".bed", sep = ""), length(clusters.q.list))
		}else{
			track.file <- paste(sample.labels, ".tagClusters", r, ".bed", sep = "")
		}
		
	}else if(what == "consensusClusters"){
		
		clusters <- object@consensusClusters
		colnames(clusters)[1] = "cluster"
		if(!(colorByExpressionProfile)){
			cols <- as.list(rep("0,0,0", length(sample.labels)))
		}else{
			clustering.method <- object@consensusClustersExpressionClusteringMethod		
		cl <- object@consensusClustersExpressionClasses
		n <- names(cl)
		cl <- .extract.cluster.info(cl = cl)
		cl <- data.frame(cluster = n, x.cor = cl[,1], y.cor = cl[,2])		
		clusters <- merge(clusters, cl, by.x = "cluster", by.y = "cluster")
		}
		
		if(length(qLow) > 0 & (length(object@consensusClustersQuantileLow)>0 & length(object@consensusClustersQuantileUp)>0)) {

		if(paste("q_", qLow, sep = "") %in% colnames(object@consensusClustersQuantileLow[[1]]) & paste("q_", qUp, sep = "") %in% colnames(object@consensusClustersQuantileUp[[1]])){
		
		use.blocks <- T
		q.low <- object@consensusClustersQuantileLow
		q.low <- lapply(q.low, function(x) {colnames(x)[2:ncol(x)] <- paste("qLow_", do.call(rbind, strsplit(colnames(x)[2:ncol(x)], split = "_", fixed = T))[,2], sep = ""); return(x)})
		q.up <- object@consensusClustersQuantileUp
		q.up <- lapply(q.up, function(x) {colnames(x)[2:ncol(x)] <- paste("qUp_", do.call(rbind, strsplit(colnames(x)[2:ncol(x)], split = "_", fixed = T))[,2], sep = ""); return(x)})
		q <- lapply(as.list(1:length(sample.labels)), function(x) {merge(q.low[[x]], q.up[[x]], by.x = "cluster", by.y = "cluster")})
		names(q) <- sample.labels
						
		cumsums <- object@CTSScumulativesConsensusClusters
		dom.pos <- list()
		for(i in 1:length(cumsums)){
			a <- lapply(cumsums[[i]], function(y) {.get.dominant.ctss(as.numeric(y), isCumulative = T)})
			b <- data.frame(cluster = as.integer(names(a)), dominant_ctss = unlist(a))
			dom.pos[[i]] <- b
		}
		names(dom.pos) <- names(cumsums)
		
		clusters.q.list <- lapply(as.list(sample.labels), function(x) {a <- merge(clusters, dom.pos[[x]], by.x = "cluster", by.y = "cluster", all.x = F, all.y = T); a$dominant_ctss <- a$start + a$dominant_ctss; b <- merge(a, q[[x]], by.x = "cluster", by.y = "cluster"); return(b)})
		track.names <- paste(sample.labels, paste(" (consensus clusters q(", qLow, ")-q(",qUp,"))", sep = ""), sep = "")
		r <- paste(".qLow", qLow, "_qUp", qUp, sep = "")
		quantiles <- T
		}else{
			stop("No data for given quantile positions! Run 'quantilePositions()' function for desired quantiles first, or omit 'qLow' and 'qUp' parameters to use start and end coordinates instead!")
		}
		}else if(length(qLow) > 0){
			stop("No data for given quantile positions! Run 'quantilePositions()' function for desired quantiles first, or omit 'qLow' and 'qUp' parameters to use start and end coordinates instead!")
		}else{
			use.blocks <- F
			clusters.q.list <- list(clusters)			
			track.names <- "Consensus clusters"
			r <- ""
			oneFile <- T
			quantiles <- F
		}
				   
		if(oneFile){
			if(quantiles){
			track.file <- rep(paste("All.samples.consensusClusters", r, ".bed", sep = ""), length(clusters.q.list))
			}else{
				track.file <- "ConsensusClusters.bed"
			}
		}else{
			track.file <- paste(sample.labels, ".consensusClusters.", r, ".bed", sep = "")		
		}
		
	}else{
		stop("'what' parameter must be one of the (\"CTSS\", \"tagClusters\", \"consensusClusters\")")
	}

	if(colorByExpressionProfile == TRUE){
		
		itemRgb <- TRUE
		cols.init <- c("red", "gold", "green", "blue")
		color.matrix.solid <- .myColorMatrix(matrix(cols.init, nrow = 2), nrows = max(cl$y.cor)+1, ncols = max(cl$x.cor)+1)
		color.matrix.solid <- matrix(as.vector(color.matrix.solid), ncol = max(cl$x.cor)+1, byrow = F)
		cols <- lapply(clusters.q.list, function(x) {m <- as.matrix(x[,c("y.cor", "x.cor")]); a <- apply(t(col2rgb(color.matrix.solid[m+1])), 1, function(x) {paste(x, collapse = ',')}); return(a)})
		if(clustering.method == "som"){
			names <- lapply(clusters.q.list, function(x) {paste(x[,"x.cor"], x[,"y.cor"], sep = "_")})
		}else if(clustering.method == "kmeans"){
			names <- lapply(clusters.q.list, function(x) {x[,"x.cor"]})
		}
		
	}else{
		itemRgb <- FALSE
		names <- as.list(rep(".", length(clusters.q.list)))
	}
	
	track.descriptions <- track.names
	
	for(i in 1:length(clusters.q.list)){
		if(i == 1){
			app <- F
		}else{
			app <- T
		}
		if(!oneFile){
			if(file.exists(track.file[i])){
				file.remove(track.file[i])
			}				
		}
		.make.cluster.bed.track(clusters.q = clusters.q.list[[i]], use.blocks = use.blocks, q.low = qLow, q.up = qUp, track.file = track.file[i], track.name = track.names[i], track.description = track.descriptions[i], cols = cols[[i]], name = names[[i]], itemRgb = itemRgb, app = app)
		
	}
	
	if(oneFile){
		message("\nFile '", track.file[1], "' has been created in your working directory (", getwd(), ")")	
	}else{
		message("\nFiles '", sub(sample.labels[1], "*", track.file[1]), "' for all samples have been created in your working directory (", getwd(), ")")
	}
	
}
)





