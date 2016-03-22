######################################################################################
# mvBM function, modified from the code provided in 'evomap' (Smaers et al., 2016).
# We emphasize that this code produces identical results to the mvBM function from 
# Jeroen Smaers' 'evomap' package for a given sigma2.
#
# There are five modifications to the original code:
#
# First, the initial estimation of sigma2 is done within the function using maximum likelihood, 
# such that the user doesn't have to do this in a separate step and feed the output to the mvBM.  
#
# Second, the output of the function is now the transformed phylogeny rather than just
# branch lengths, such that the user does not need to set the branch lengths of their
# tree equal to the mvBM output in a separate step. 
#
# Third, quite a bit of unnecessary and redundant code has been removed to enhance
# readability of the script. There was some whacky stuff in the original code, and this
# modified code performs more functions with more annotation in fewer lines.
#
# Fourth, a consistent assignment operator is used throughout the code (<-). This is a 
# purely aesthetic change, because I had to look at this code for a long time and I
# couldn't handle the random switching back and forth between assignment operators. I 
# honestly don't care whether people use <- or =, but for heavens sake pick ONE.
#
# Finally, the code has been heavily annotated to clarify what is going on in various
# steps.
######################################################################################

mvBM <- function(data, tree) {
	
  # Get sigma2 using restricted ML
	sigma2 <- ace(data, tree, method="REML")
	sigma2 <- sigma2$sigma2[1]
 		
  # Make phylo matrix and list nodes
  phy.matrix <- data.frame(tree$edge, tree$edge.length, data[tree$edge[,2]], row.names=(1:nrow(tree$edge)))
  names(phy.matrix) <- c("Anc","Desc","Length","Value")
  nodes_extant <- 1:length(data)
  nodes_extinct <- (length(data)+1):(2*(length(data))-1)

  # Compute "global" estimate
  Pk <- c()

	  # Loop through internal nodes
      for(j in nodes_extinct) {
      	# Loop through terminal nodes, compute nominator and denominator of Pk
        nominator <- denominator <- c() 
        for(i in nodes_extant){
             nominator <- rbind(nominator,(data[i]/(dist.nodes(tree)[i,j]^2)))
             denominator <- rbind(denominator,(1/(dist.nodes(tree)[i,j]^2)))
         }
      	 # Save a Pk for each internal node
         Pk <- rbind(Pk, sum(nominator)/sum(denominator))
       }

	  # Make data frame w Pk values
      Pk <- data.frame(value=Pk, nodes=nodes_extinct, row.names=1:length(nodes_extinct))

	  # Rescaled branch lengths
      rBL <- rBL_edge <- c()

	  # Loop through internal nodes in reverse order (traverse tree from tips to root)
    for(i in rev(nodes_extinct)) {
      # Get indices of i's descendent branches in phy.matrix
      sister_branches <- which(phy.matrix$Anc==i)
		  # Get trait value of i's descendents
      X1 <- phy.matrix$Value[sister_branches[1]]
      X2 <- phy.matrix$Value[sister_branches[2]]
		  # Get Pk of i
      Pk_anc <- Pk$value[which(Pk$nodes==i)]
		  # Take the average of Pk and the trait values of i's descendents
      Ax <- (Pk_anc+X1+X2)/3
		  # Compute new branch lengths
      T1 <- ((X1-Ax)^2)/sigma2 + phy.matrix$Length[sister_branches[1]]
      T2 <- ((X2-Ax)^2)/sigma2 + phy.matrix$Length[sister_branches[2]]
		  # Store estimated trait value for internal node i in phy.matrix
      phy.matrix$Value[which(phy.matrix$Desc==i)] <- Ax
		  # Record new branch lengths
      rBL <- c(rBL, T1, T2)
      rBL_edge <- c(rBL_edge, sister_branches)
    }
       
  # Return tree
  names(rBL) <- rBL_edge
  tree$edge.length <- rBL[order(as.numeric(names(rBL)))]
  return(tree)
}
