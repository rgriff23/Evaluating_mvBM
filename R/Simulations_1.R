# This code simulates a three taxa tree with a burst on the longest branch. The focus
# is on how ML, PIC, and mvBM reconstruct the one non-root internal node.
# The expectation is that ML will spread change across the branches such that the
# internal node is estimate closer to the tip value at the end of the long branch,
# while PIC will be completely insensitive to the burst on the long branch.
# mvBM will likely be intermediate, since it is not entirely blind to the burst branch,
# but is biased towards a local ML reconstruction.
##########################################################################################
# PREPARATIONS
##########################################################################################

# Load necessary packages and functions
library(phytools)
source('mvBM.R') # user must firt navigate to directory w all the R files

##########################################################################################
# SIMULATIONS
##########################################################################################

# Set seed
set.seed(23)

# Create tree
N <- 3
tree <- tree.b <- pbtree (n = N)

# Define burstbranch and internal node of interest
burstbranch <- 4
node <- 5

# Simulate burst 100X on burstbranch
nsim <- 200
sigma2 <- 0.01
sims <- fastBM(tree.b, sig2=sigma2, nsim=nsim, internal=TRUE)
burst <- seq(0.01, 1.5, length.out=nsim)
sims[3,] <- burst

# Separate simulated tip values (dat) and node values (anc)
dat <- sims[1:N,]
anc <- sims[as.character(node),]

# Amount of simulated change on burst branch
burst <- sims["t1",]

##########################################################################################
# ANCESTRAL STATE RECONSTRUCTIONS
##########################################################################################

# Store internal node estimate
ml5 <- pic5 <- mvBM5 <- c()

# Loop through sims and do reconstructions
for (i in 1:nsim) {
  # Maximum likelihood
  ml <- ace(dat[,i], tree, method="REML")
  ml5[i] <- ml$ace[2]
  # PIC
  pic <- ace(dat[,i], tree, method="pic")
  pic5[i] <- pic$ace[2]
  # mvBM
  mvtree <- mvBM(dat[,i], tree)
  mv <- ace(dat[,i], mvtree, method="REML")
  mvBM5[i] <- mv$ace[2]  
}

##########################################################################################
# PLOTS
##########################################################################################

# Inaccuracy of reconstructions
ml.err <- ml5 - anc
pic.err <- pic5 - anc
mv.err <- mvBM5 - anc

# Phylogeny and line plot to show inaccuracy as burst size increases
layout(matrix(1:2, 1, 2))
plot(tree)
nodelabels(frame="circle", bg="white", font=2)
edgelabels(text=c("1 ", "2 ", "3 ", "4 "), frame="rect", bg="black", col="white", font=2)
mtext("A", adj=0, cex=2, line=1)
plot(ml.err ~ burst, type="l", ylim=c(-0.1, 0.5), ylab="Error of estimate at node 5", xlab="Size of burst on branch 4", col="red")
lines(pic.err ~ burst, col="blue")
lines(mv.err ~ burst, col="purple")
mtext("B", adj=0, cex=2, line=1)
legend(1.15, 0.675, legend=c("cvBM", "mvBM", "PIC"), col=c("red", "purple", "blue"), lwd=3, xpd=TRUE, cex=0.6)

##########################################################################################
# END
##########################################################################################