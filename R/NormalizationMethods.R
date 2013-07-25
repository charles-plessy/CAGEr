######################################################################################
# Funtions for normalizing CAGE tag count to a referent power-law distribution
# Reference: Balwierz, P. J., Carninci, P., Daub, C. O., Kawai, J., Hayashizaki, 
# Y., Van Belle, W., Beisel, C., et al. (2009). Methods for analyzing deep sequencing 
# expression data: constructing the human and mouse promoterome with deepCAGE data. 
# Genome Biology, 10(7), R79.
#


setGeneric(
name="normalizeTagCount",
def=function(object, method = "powerLaw", fitInRange = c(10, 1000), alpha = 1.25, T = 10^6){
	standardGeneric("normalizeTagCount")
}
)

setMethod("normalizeTagCount",
signature(object = "CAGEset"),
function (object, method, fitInRange, alpha, T){
	
	objName <- deparse(substitute(object))
	
	message("\nNormalizing tag count...")
	sample.labels <- sampleLabels(object)
	tag.count <- object@tagCountMatrix

	if(method == "powerLaw"){
	first <- TRUE

	for(x in sample.labels) {
		
		message("\t-> ", x)
		fit.coef <- .fit.power.law.to.reverse.cumulative(values = as.integer(tag.count[,x]), val.range = fitInRange)
		norm.tag.count <- .normalize.to.reference.power.law.distribution(values = tag.count[,x,drop=F], lin.reg.coef = fit.coef, alpha = alpha, T = T)
		if(first == TRUE) {
			ctss.all.norm <- data.frame(norm.tpm = norm.tag.count)
		}else{
			ctss.all.norm <- cbind(ctss.all.norm, norm.tag.count)
		}
		
		first <- FALSE
		
	}
	
	}else if(method == "simpleTpm"){
		ctss.all.norm <- as.data.frame(apply(tag.count, 2, function(x) {x/sum(x) * 10^6}))
	}else if(method == "none"){
		ctss.all.norm <- as.data.frame(tag.count)
	}else{
		stop("'method' must be one of the (\"powerLaw\", \"simpleTpm\", \"none\")")
	}
	
	colnames(ctss.all.norm) = sample.labels
	
	object@normalizedTpmMatrix <- ctss.all.norm
	assign(objName, object, envir = parent.frame())
	invisible(1)
	
}
)
