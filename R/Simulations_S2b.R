# This code shows that likelihoods improve when mvBM is looped multiple times
##########################################################################################
# LOAD PACKAGES AND FUNCTIONS
##########################################################################################

# Load necessary packages and functions
library(phytools)
library(geiger)
source('mvBM.R') # user must firt navigate to directory w all the R files

##########################################################################################
# SIMULATE TREES AND DATA
##########################################################################################

# Set seed
set.seed(23)

# Generate pure birth tree with 100 tips
ntips <- 100
tree <- pbtree(n=ntips)

# Simulate 1 continuous BM trait with sigma2=0.01
sigma2 <- 0.01
bm.tips <- fastBM(tree, sig2=sigma2, nsim=1)

##########################################################################################
# FIT cvBM (regular tree), mvBM (mvBM tree), and ML-mvBM
##########################################################################################

lk <- c()
temp <- tree
loops <- 10
for (i in 1:loops) {
  temp <- mvBM(bm.tips, temp)
  lk[i] <- logLik(fitContinuous(temp, bm.tips))
  print(i)
}

##########################################################################################
# PLOTS
##########################################################################################

# Compare likelihoods with line plot (Figure S2b)
quartz()
plot(lk ~ c(1:loops), type="l", col="black", xlab="Iteration", ylab="Log Likelihood")

##########################################################################################
# END
##########################################################################################
