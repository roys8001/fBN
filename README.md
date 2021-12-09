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
