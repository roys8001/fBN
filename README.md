
---
output:
  pdf_document: default
  html_document: default
---
## Functional Bayesian Network - Intermediate Project Update

# Introduction

Continuously monitoring multiple interrelated variables is routinely conducted in longitudinal studies, constituting a type of multivariate functional data. One critical task is to explore the dependence/ causal structure among the targets of interests. In this package we tend to develop a Bayesian approach to learn semi-parametrically directed acyclic graphical(DAG) models for multivariate functional data. Specifically we allow for deviation from Gaussian process of functions which facilitates unique structure identification and robustness against outliers. The model builds on a dirichlet process scale mixtures in the graphical modeling of principal scores and a data driven approach in learning principal functions which allows for flexibility in capturing within and between function function covariance. The posterior inference is carried out with Markov Chain Monte Carlo. 

# Installation instructions

Run the following command - devtools::install_github("roys8001/fBN")

# To Do

The following things are left to be done:

1) Build R wrapper.

2) Add three functions to find principal scores, sampling  error and to generate causal structure information to the fBN.cpp file.

3) Add a function to obtain the initial parameters to the fBN.R file.

4) Perform compatibility checks on user supplied input.

5) If time permits, create a vignette to demonstrate the usage.


