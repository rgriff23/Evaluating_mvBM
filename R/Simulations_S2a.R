# This code shows that mvBM trees fit the data better than the original tree even when
# data are simulated under cvBM.
##########################################################################################
# LOAD PACKAGES AND FUNCTIONS
##########################################################################################

# Load necessary packages and functions
library(phytools)
library(geiger)
source('mvBM.R') # user must firt navigate to directory w all the R files
source('ML_mvBM.R', chdir = TRUE)

##########################################################################################
# SIMULATE TREES AND DATA
##########################################################################################

# Set seed
set.seed(23)

# Generate pure birth tree with 100 tips
ntips <- 100
tree <- pbtree(n=ntips)
#quartz(); plot(tree); nodelabels(); edgelabels()

# Simulate 100 continuous BM traits with sigma2=0.01
nsims <- 100
sigma2 <- 0.01
bm.tips <- fastBM(tree, sig2=sigma2, nsim=nsims)

##########################################################################################
# FIT cvBM (regular tree), mvBM (mvBM tree), and ML-mvBM
##########################################################################################

lk.cvBM <- lk.mvBM <- lk.MLmvBM <- c()

for (i in 1:nsims) {
  # Fit original tree
  lk.cvBM[i] <- logLik(fitContinuous(tree, bm.tips[,i]))
  # Fit mvBM to BM data
  mvtree <- mvBM(bm.tips[,i], tree)
  lk.mvBM[i] <- logLik(fitContinuous(mvtree, bm.tips[,i]))
  # Fit ML_mvBM to BM data
  MLmvtree <- ML_mvBM(bm.tips[,i], tree)
  lk.MLmvBM[i] <- logLik(fitContinuous(MLmvtree, bm.tips[,i]))
  # Monitor progress
  print(i)
}

##########################################################################################
# PLOTS
##########################################################################################

# Compare likelihoods with boxplot (Figure S2a)
quartz()
models <- c(rep("cvBM", length(lk.cvBM)), rep("mvBM", length(lk.cvBM)), rep("ML-mvBM", length(lk.cvBM)))
models <- factor(models, levels=c("cvBM", "mvBM", "ML-mvBM"))
d <- data.frame(loglik=c(lk.cvBM, lk.mvBM, lk.MLmvBM), models=models)
boxplot(loglik ~ models, data=d, col="gray", ylab="Log likelihood", xlab="Model")

##########################################################################################
# END
##########################################################################################
