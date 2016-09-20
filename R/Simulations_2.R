# This code demonstrates error in ancestral state estimates near the site of a burst.
# Run time is 2.5-3 hrs on my Macbook Pro (2.6 GHz Intel Core i5, 8 GB 1600 MHz DDR3)
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

# Generate a pure birth tree with 100 tips
ntips <- 100
tree <- pbtree(n=ntips)
quartz(); plot(tree); nodelabels(); edgelabels()

# Simulate 5 continuous BM traits with sigma2=0.01 for each branch
nsims <- 5
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

# Median range of simulated values?
spread <- c()
for (i in 1:(length(tree$edge.length)*nsims)) {
  spread[i] <- abs(diff(range(burst.sims[,i])))
}
median(spread)

# Record trait change and error for node at base of each burst branch
burst <- c()
err <- c()
for (i in 1:(length(tree$edge.length)*nsims)) {
  burst[i] <- diff(burst.sims[tree$edge[burstbranch[i],],i])
  actual <- burst.anc[as.character(tree$edge[burstbranch[i],1]),i]
  estimated <- burst.mvBM[[i]]$ace[as.character(tree$edge[burstbranch[i],1])]
  err[i] <- diff(c(actual, estimated))
}
range(err)

# Plot error as a function of burst size
plot(err ~ burst, ylab="Error at base of burst branch", xlab="Trait change on burst branch")

##########################################################################################
# END
##########################################################################################
