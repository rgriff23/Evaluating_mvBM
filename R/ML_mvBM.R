# Return appropriately transformed phylogeny
ML_mvBM <- function (data, tree) {
	
	# Reconstruct ancestral states using restricted ML
	ace.out <- ace(data, tree, method="REML")
	
	# Compute directional contrasts and transform tree
	nodes <- c(data, ace.out$ace)
	names(nodes) <- 1:length(nodes)
	tree$edge.length <- abs(nodes[tree$edge[,2]] - nodes[tree$edge[,1]])
	
	# Return tree
	return(tree)
}