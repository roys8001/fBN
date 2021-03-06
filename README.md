Functional Bayesian Network - fBN
================
SAPTARSHI ROY
05/12/2021

## Introduction

Continuously monitoring multiple interrelated variables is routinely
conducted in longitudinal studies, constituting a type of multivariate
functional data. One critical task is to explore the dependence/ causal
structure among the targets of interests. In this package we tend to
develop a Bayesian approach to learn semi-parametrically directed
acyclic graphical(DAG) models for multivariate functional data.
Specifically we allow for deviation from Gaussian process of functions
which facilitates unique structure identification and robustness against
outliers. The model builds on a dirichlet process scale mixtures in the
graphical modeling of principal scores and a data driven approach in
learning principal functions which allows for flexibility in capturing
within and between function function covariance. The posterior inference
is carried out with Markov Chain Monte Carlo.

## Installation instructions

``` r
#devtools::install_github("roys8001/fBN")
```

## Example

``` r
################################################################################
########################### generate simulation data ###########################
################################################################################

# library(rmutil); library(orthogonalsplinebasis); library(igraph)
# 
# numx = 500; numc = 10; numz = 2; numk = 4; numt = 200; numxt = 10; vstrength = 1; vnoise = 0.5
# 
# # simulation truth of outer graph structure: DAG
# G0 = matrix(0, numc, numc); G0[lower.tri(G0)] = as.numeric(runif(sum(lower.tri(G0))) < 2 / numc)
# # simulation truth of inner graph structure
# GS = matrix(0, numc * numk, numc * numk)
# for(j in 2 : numc) {
#   for(l in 1 : (j - 1)) {
#     if(G0[j, l] == 1) {
#       cindex = ((j - 1) * numk + 1) : (j * numk); findex = ((l - 1) * numk + 1) : (l * numk)
#       # selection of edges
#       GS[cindex, findex] = as.numeric(runif(numk * numk) < 0.5)
#       while(sum(GS[cindex, findex]) == 0) GS[cindex, findex] = as.numeric(runif(numk * numk) < 0.5)
#     }
#   }
# }
# 
# # nonzero signal strength
# BS = matrix(0, numc * numk, numc * numk)
# BS[GS != 0] = runif(sum(GS != 0), vstrength / 2, vstrength) # uniform
# BS[GS != 0] = BS[GS != 0] * (1 - 2 * (runif(sum(GS != 0)) < 0.5))
# 
# # generate principal scores
# score = matrix(rlaplace(numx * numc * numk, s = vnoise), numx, numc * numk)
# score = t(tcrossprod(solve(diag(numc * numk) - BS), score))
# 
# # generate all available observation time
# t0 = sort(runif(numt)); tx = array(0, dim = c(numx, numc, numxt))
# # generate observation time for each sample and node
# for(i in 1 : numx) for(j in 1 : numc) tx[i, j, ] = sample(t0, numxt)
# # summarize observation time for each node
# tobs = NULL; for(j in 1 : numc) tobs = c(tobs, sort(unique(as.vector(tx[, j, ]))))
# 
# # generate functional principal components
# lb = min(t0); ub = max(t0); ib = quantile(t0, c(0, 1))
# fPC = evaluate(OBasis(knots = c(rep(lb, numk - 1), ib, rep(ub, numk - 1)), order = numk), t0)
# 
# # generate baseline variables
# zobs = matrix(rnorm(numx * numz), numx, numz)
# # polynomial basis for baseline effects
# zbase = bs(t0); nums = ncol(zbase); ecoef = matrix(runif(numc * numz * nums, vstrength / 2, vstrength), numc, numz * nums)
# ecoef = ecoef * (1 - 2 * (runif(length(ecoef)) < 0.5))
# 
# # generate observation
# xobs = matrix(0, numx, length(tobs)); cindex = c(1, which(tobs[-1] < tobs[-length(tobs)]) + 1, length(tobs) + 1); e0 = rep(0, numc)
# for(j in 1 : numc) {
#   # node index
#   fposition = ((j - 1) * numk + 1) : (j * numk); cposition = cindex[j] : (cindex[j + 1] - 1)
#   # calculate signal strength
#   ee = tcrossprod(score[, fposition], fPC[match(tobs[cposition], t0), ])
#   for(l in 1 : numz) ee = ee + crossprod(t(zobs[, l]), tcrossprod(ecoef[j, ((l - 1) * nums + 1) : (l * nums)], zbase[match(tobs[cposition], t0), ]))
#   e0[j] = vnoise * mean(abs(ee))
#   # match sample with observation time
#   for(i in 1 : numx) {
#     xposition = match(tx[i, j, ], tobs[cposition]); tposition = match(tx[i, j, ], t0)
#     # generate observed data
#     xobs[i, cposition[xposition]] = e0[j] * rnorm(numxt) + crossprod(ecoef[j, ], kronecker(zobs[i, ], t(zbase[tposition, ]))) + tcrossprod(score[i, fposition], fPC[tposition, ])
#   }
# }
# # fill zeros with NA
# xobs[xobs == 0] = NaN
# 
# ################################################################################
# ################################ MCMC inference ################################
# ################################################################################
# 
# source('fBN.R')
# Rcpp::sourceCpp('fBN.cpp')
# 
# init_para = get_initial_parameter(xobs, zobs, tobs)
# # index information
# cindex = init_para$cindex; findex = init_para$findex; tindex = init_para$tindex
# # spline information
# sbasis = init_para$sbasis; penat = init_para$penat
# # parameter information
# fscore = init_para$fscore; scoef = init_para$scoef; lambda = init_para$lambda
# zcoef = init_para$zcoef; eta = init_para$eta; svar = init_para$svar; numk = init_para$numk
# # causal information
# gstructure = init_para$gstructure; fstructure = init_para$fstructure; mass = init_para$mass
# fcoef = init_para$fcoef; fvar = init_para$fvar; cvar = init_para$cvar; eprob = init_para$eprob
# cmat = init_para$cmat; degree = init_para$degree
# 
# RECORDZ = NULL; RECORDS = NULL; RECORDF = NULL; RECORDV = NULL; RECORDG = NULL
# 
# iter = 1; maxit = 1000; burnin = 500; interval = 10
# while(iter <= maxit) {
#   # baseline effects
#   zinfo = get_zcoef(xobs, zobs, scoef, sbasis, fscore, eta, svar, tindex, cindex, findex)
#   zcoef = zinfo$zcoef; eta = zinfo$eta
#   # spline information
#   sinfo = get_scoef(xobs, zobs, zcoef, scoef, sbasis, penat, fscore, svar, lambda, tindex, cindex, findex)
#   scoef = sinfo$scoef; lambda = sinfo$lambda
#   # principal scores
#   fscore = get_fscore(xobs, zobs, zcoef, scoef, sbasis, fvar, fcoef, svar, tindex, cindex, findex)
#   # sampling error
#   svar = get_svar(xobs, zobs, zcoef, scoef, sbasis, fscore, tindex, cindex, findex)
#   # causal structure information
#   ginfo = get_ginfo(fscore, fstructure, gstructure, fvar, cmat, degree, mass, findex, cvar, eprob)
#   fvar = ginfo$fvar; fcoef = ginfo$fcoef; gstructure = ginfo$gstructure; fstructure = ginfo$fstructure
#   cvar = ginfo$cvar; eprob = ginfo$eprob; cmat = ginfo$cmat; degree = ginfo$degree; mass = ginfo$mass
#   
#   if(iter %in% seq(burnin, maxit, interval)) {
#     RECORDZ = c(RECORDZ, list(zinfo))
#     RECORDS = c(RECORDS, list(sinfo))
#     RECORDF = c(RECORDF, list(fscore))
#     RECORDV = c(RECORDV, list(svar))
#     RECORDG = c(RECORDG, list(ginfo))
#   }
#   
#   iter = iter + 1
# }
# 
# save(RECORDZ, RECORDS, RECORDF, RECORDV, RECORDG, file = "RECORD.RData")
```

# Summarizing the results

``` r
################################################################################
########################## summarize simulation result #########################
################################################################################

# posterior estimate of graph structure
# GEST = matrix(0, numc, numc)
# for(i in 1 : length(seq(burnin, maxit, by = interval))) GEST = GEST + RECORDG[[i]]$gstructure
# GEST = matrix(as.numeric(GEST / length(seq(burnin, maxit, by = interval)) > 0.5), ncol = numc)
# 
# par(mfrow = c(1, 2)); par(mar = c(2, 2, 2, 2))
# 
# plot(graph_from_adjacency_matrix(t(G0)), main = "simulation truth")
# plot(graph_from_adjacency_matrix(t(GEST)), main = "posterior estimate")
# 
# tpos = sum(GEST[G0 == 1]); tneg = sum(1 - GEST[G0 == 0])
# fpos = sum(GEST[G0 == 0]); fneg = sum(1 - GEST[G0 == 1])
# # accuracy summary true positive / false discovery
# TPR = tpos / (tpos + fneg); FDR = fpos / (fpos + tpos)
# # Matthews correlation coefficient
# MCC = (tpos * tneg - fpos * fneg) / sqrt((tpos + fpos) * (tpos + fneg) * (tneg + fpos) * (tneg + fneg))
# 
# par(mfrow = c(3, 4)); par(mar = c(2, 2, 2, 2))
# # estimated PC components
# ePC = crossprod(t(sbasis), scoef)
# for(j in 1 : ncol(fPC)) plot(fPC[, j], type = "l", xlab = "", ylab ="", main = "fPC")
# for(j in 1 : ncol(ePC)) plot(ePC[, j], type = "l", xlab = "", ylab ="", main = "ePC")
# 
# for(j in 1 : numc) {
#   for(l in 1 : numz) {
#     tbaseline = as.vector(crossprod(ecoef[j, ((l - 1) * nums + 1) : (l * nums)], t(zbase)))
#     ebaseline = as.vector(crossprod(zcoef[j, ((l - 1) * ncol(sbasis) + 1) : (l * ncol(sbasis))], t(sbasis)))
#     # true baseline coefficient
#     plot(tbaseline, type = "l", xlab = "", ylab ="", main = "truth")
#     # estimated baseline coefficient
#     plot(ebaseline, type = "l", xlab = "", ylab ="", main = "estimate")
#   }
# }
```
