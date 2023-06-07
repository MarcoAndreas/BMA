
# Being Bayesian about 
# Learning Gaussian Bayesian networks
# from incomplete data

# R code for the Bayesian Model Averaging (BMA) approach
# for Gaussian Bayesian networks with missing data

# paper to appear in the
# International Journal of Approximate Reasoning

# author: Marco Grzegorczyk (Groningen University, NL)
# Email to: m.a.grzegorczyk@rug.nl

#########################################################################

# Two libraries are needed to run the code
# Both libraries are available from CRAN (and must be installed in advance)
library(pracma)
library(MASS)

# Specify your working directory (where the two files with R code files have been saved)
# MCMC_ALGORITHMS.R
# MORE_FUNCTIONS.R

# For example
setwd("X:/BMA")

# Read the two provided files with R codes for Bayesian Model Averaging (BMA)
source("MCMC_ALGORITHM.R")
source("MORE_FUNCTIONS.R")

#################################################################

# For illustration purposes:
# Generate a small synthetic data set 
# with n=5   nodes 
# and  N=100 observations

n = 5
N = 100

data = matrix(0,n,N)

# Just Gaussian random numbers
data[1,] = rnorm(N)
data[2,] = rnorm(N)

# Node X3 depends on nodes X1 and X2 (X3 <- X1,X2)
data[3,] = data[1,] + data[2,] + 0.1 * rnorm(N)

# Just Gaussian random numbers
data[4,] = rnorm(N)

# Node X5 depends on node X4 (X5 <- X4)
data[5,] = data[4,] + 0.1 * rnorm(N)

#################################################################

# Delete a fraction of 'p_delete' data points

p_delete = 0.1  

n_datapoints = prod(dim(data))
  
n_delete     = ceiling(n_datapoints * p_delete)
  
ind = sample(1:n_datapoints,n_delete,replace=FALSE)

data[ind] = NA

##################################################################

# Optionally load your own incomplete data set 'data'
# It must be a data matrix, where 
# data[i,j] is the j-th value of node Xi

# And the data are supposed to have missing ('NA') values 

##################################################################

# Specify an initial DAG (here a DAG without any edges)

DAG = matrix(0,n,n)

# Total number of MCMC iterations
T_mcmc = 10000

# Specify a thin-out factor
f_thin = 100

# No. of samples will be T_mcmc * (1/f_thin) (here: 100)

# Specify MCMC tuning parameter 
# q is the probability that missing data will be re-sampled
q = 0.05

# Run the MCMC simulation
out = strMCMC_NOVEL(data,DAG0,T_mcmc,f_thin,q)
  
# out[[1]][[1]] is the adjacency matrix of the initial DAG
# out[[2]][[1]] is the BGe score        of the initial DAG

# out[[1]][[k+1]] is the adjacency matrix of the k-th sampled DAG
# out[[2]][[k+1]] is the BGe score        of the k-th sampled DAG

# Process the results by computing CPDAGs and averaging across them
# Skip the first 50 DAGs (to take burn in phase into account)

results = cpdag_list(out,50)
  
######################################################################

# Take out the edge scores

BMA_SCORES = results[[3]]

# Edge scores are given by:

BMA_SCORES

#######################################################################

# If the true DAG is known:

# Here the adjacency matrix of the true DAG is:

TRUE_DAG      = matrix(0,n,n)
TRUE_DAG[1,3] = 1
TRUE_DAG[2,3] = 1
TRUE_DAG[4,5] = 1

# One can compute the achieved AUROC and rSHD score

auroc_score = compute_AUROC(BMA_SCORES,TRUE_DAG)

rSHD_score  =   compute_SHD(BMA_SCORES,TRUE_DAG) 

#######################################################################

rSHD_score

auroc_score
  

