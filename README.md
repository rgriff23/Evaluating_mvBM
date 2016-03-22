# Evaluating the legitimacy of the  mvBM method by Smaers et al. (2016)

[Randi H. Griffin](http://rgriff23.github.io/) and [Gabriel S. Yapuncich](http://www.gabrielyapuncich.com/)

___

This is a project to show that Smaers et al.'s (2016) new phylogenetic comparative method is problematic. A manuscript is in prep.

## Reproducible simulation study


The **R** folder contains the data and R code needed to replicate our simulation study. Several R packages must be installed for this code to work.

```
install.packages("ape")
install.packages("geiger")
install.packages("phytools")
```

Six separate simulation studies are included in our study. These simulations and associated figures can be replicated with the following code (assuming the R folder is your working directory).

```
# Figure 1
source("Simulations1_overfitting.R")

# Figure 2
source("Simulations2_relativeaccuracy.R")

# Figure 3
source("Simulations3_errorbias.R")

# Figure 4
source("Simulations4_30taxatree.R")

# Supplement S2
source("SimulationsS1_10iterations.R")

# Supplement S3
source("SimulationsS2_5burst.R")
```

The function `mvBM.R` is our implementation of Jeroen Smaers' mvBM function from the evomap package, which can be found [here](https://github.com/rgriff23/evomap/blob/master/R/mvBM.R). The function `ML_mvBM.R` is our modified version of mvBM, which is used in the first set of simulations (`Simulations1_overfitting.R`). 

## References

- Smaers JB, Mongle CS, Kandler A. 2016. A multiple variance Brownian motion framework for estimating variable rates and inferring ancestral states. Biol J Linn Soc, DOI: 10.1111/bij.12765
 
___
