#####
# Function that calculates position of quantiles for CTSS clusters based on distribution of tpm within the cluster
# ARGUMENTS: cluster.cumsums - named list of vectors containing cumulative sum for each cluster (returned by 'get.cumsum' function)
#            coors - data frame of clusters containing least one column, *start (genomic start position of the cluster)
#			 q.low - desired lower quantile - single value or a vector of values
#            q.up - desired upper quantile - single value or a vector of values (must be the same length as q.low)
# RETURNS: 


.get.quant.pos <- function(cluster.cumsums, coors, q = NULL, useMulticore = FALSE, nrCores = NULL) {
	
  getQuantilepos <- function(quant, cumsum) length(Filter(isTRUE, cumsum/max(cumsum) < quant)) + 1
  useMulticore <- .checkMulticore(useMulticore)

	if(length(q) > 0) {
		if(useMulticore == TRUE){
			cluster.q = mclapply(cluster.cumsums, function(x) {sapply(q, getQuantilepos, x)}, mc.cores = .getNrCores(nrCores))
		}else{
			cluster.q = lapply(cluster.cumsums, function(x) {sapply(q, getQuantilepos, x)})
		}
		cluster.q = do.call(rbind, cluster.q)
		if(class(cluster.q[1,1]) == 'list') {
			un = unlist(lapply(cluster.q, function(x) {length(unlist(x))}))
			l = lapply(cluster.q, function(x) {unlist(x)})
			a1 = do.call(rbind, l[which(un !=0)])
			cluster.q = as.data.frame(matrix(rep(coors$start[which(coors$tpm>0)], ncol(cluster.q)), ncol = ncol(cluster.q), byrow = F) -1 + as.matrix(a1))
		}else{
			cluster.q = as.data.frame(matrix(rep(coors$start[which(coors$tpm>0)], ncol(cluster.q)), ncol = ncol(cluster.q), byrow = F) -1 + cluster.q)
		}
		colnames(cluster.q) = paste('q_', q, sep = '')
		cluster.q = cbind(cluster = as.integer(names(cluster.cumsums)[which(coors$tpm>0)]), cluster.q)
		cluster.q = merge(coors, cluster.q, all.x = T, all.y = F)
	}
		
	return(cluster.q[order(cluster.q$cluster),])	
	
}
