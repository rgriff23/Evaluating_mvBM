
# This code shows that the R^2 is lower when fewer taxa (30 vs. 100) are simulated
# and also that the R^2 of mvBMestimated vs. true ancestral states varies across
# the branches of the tree, with branches near the tips exhibiting lower R^2 values
##########################################################################################
# PREPARATIONS
##########################################################################################

# Load necessary packages and functions
library(phytools)
source('~/Desktop/GitHub/Evaluating_mvBM/R/mvBM.R', chdir = TRUE)

##########################################################################################
# SIMULATE TREES AND DATA
##########################################################################################

# Set seed
set.seed(23)

# Generate a pure birth tree with 30 tips
ntips <- 30
tree <- pbtree(n=ntips)
#quartz(); plot(tree); nodelabels(); edgelabels()

# Simulate 100 continuous BM traits with sigma2=0.01 for each branch
nsims <- 100
sigma2 <- 0.01

# Simulate 'one-burst' scenario on each branch, where burst branch has sigma2=1
burst.sims <- c()
burstbranch <- sort(rep(1:length(tree$edge.length), nsims))
for (i in 1:(length(tree$edge.length)*nsims)) {
  temp <- tree
  temp$edge.length[burstbranch[i]] <- temp$edge.length[burstbranch[i]]*100
  burst.sims <- cbind(burst.sims, fastBM(temp, sig2=sigma2, internal=TRUE))
}

# Separate simulated tips from simulated ancestral states
burst.tips <- burst.sims[1:ntips,]
burst.anc <- burst.sims[(ntips+1):(2*ntips-1),]

##########################################################################################
# FIT mvBM (using ace REML rather than anc.Bayes for efficiency)
##########################################################################################

# Store ancestral state reconstructions
burst.mvBM <- list()

for (i in 1:ncol(burst.tips)) {
  # Fit mvBM to burst data
  mvtree <- mvBM(burst.tips[,i], tree)
  burst.mvBM[[i]] <- ace(burst.tips[,i], mvtree, method="REML")
  # Show progress
  print(i)
}

##########################################################################################
# PLOTS
##########################################################################################

# Compute R^2 for each estimated vs simulated
R2 <- c()
for (i in 1:length(burst.mvBM)) {
	r2 <- summary(lm(burst.mvBM[[i]]$ace ~ burst.anc[,i]))
	R2[i] <- r2$r.squared
	print(i)
}

# Compute median R^2 for each branch
medR2 <- c()
for (i in 1:length(tree$edge.length)) {
	medR2[i] <- median(R2[burstbranch==i])
}

# Plot distribution of R^2 and tree with R^2 for each branches
quartz()
layout(matrix(1:2, 1, 2))
hist(R2, col="gray", main="", xlab="R-squared")
abline(v=median(R2), lwd=2, lty=2)
par(mar=c(1,1,1,1))
plot(tree, label.offset=0.1)
edgelabels(round(medR2, 2), bg="black", col="white", font=2, cex=0.5)

##########################################################################################
# END
##########################################################################################
