# This code attempts to replicate results from Figs. 3-4 of Smaers et al., 2016
##########################################################################################
# PREPARATIONS
##########################################################################################

# Load necessary packages and functions
library(phytools)
source('mvBM.R') # user must firt navigate to directory w all the R files

##########################################################################################
# SIMULATE TREES AND DATA
##########################################################################################

# Set seed
set.seed(23)

# Generate a pure birth trees with 30 tips
ntips <- 30
tree <- pbtree(n=ntips)

# Simulate 1000 continuous BM traits with sigma2=0.01
nsims <- 1000
sigma2 <- 0.01

# For each tree, simulate '5-burst' scenario where 5 random branches have sigma2=1
burst.sims <- c()
for (i in 1:nsims) {
  temp <- tree
  burstbranch <- sample(1:(2*ntips - 2), size=5, replace=FALSE)
  temp$edge.length[burstbranch[i]] <- temp$edge.length[burstbranch[i]]*100
  burst.sims <- cbind(burst.sims, fastBM(temp, sig2=sigma2, internal=TRUE))
}

# Separate simulated tips from simulated ancestral states
burst.tips <- burst.sims[1:ntips,]
burst.anc <- burst.sims[(ntips+1):(2*ntips-1),]

##########################################################################################
# FIT cvBM AND mvBM (using ace REML rather than anc.Bayes for efficiency)
##########################################################################################

# Store ancestral state reconstructions
burst.mvBM <- list()
for (i in 1:ncol(burst.tips)) {
  mvtree <- mvBM(burst.tips[,i], tree)
  burst.mvBM[[i]] <- ace(burst.tips[,i], mvtree, method="REML")
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

# Plot distribution (Figure S3)
hist(R2, col="gray", xlab="R-squared", main="")
abline(v=median(R2), col="red", lwd=2)

##########################################################################################
# END
##########################################################################################

