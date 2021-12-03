### Loading required packages

library(fda)
library(splines2)
library(RGeode)

### Function to generate spline information

#' Title
#'
#' @param time
#' @param nknots
#'
#' @return
#' @export
#'
#' @examples
get_sinfo = function(time, nknots) {
  # get boundary and interior knots
  utime = sort(unique(time)); lknot = min(utime); uknot = max(utime)
  iknot = quantile(utime, seq(0, 1, length = (nknots + 2))[- c(1, (nknots + 2))])

  # create function objects for penalty calculation
  bsbasis = create.bspline.basis(c(lknot, uknot), breaks = c(lknot, iknot, uknot))
  mobasis = create.monomial.basis(c(lknot, uknot), nbasis = 2)

  # generate basis and penalty matrix
  sbasis = bs(utime, knots = iknot, degree = 3, Boundary.knots = c(lknot, uknot), intercept = TRUE)
  pbasis = eval.penalty(bsbasis, Lfdobj = 2, rng = c(lknot, uknot))
  ibasis = eval.penalty(bsbasis, Lfdobj = 0, rng = c(lknot, uknot))

  # transform to align with the linear/nonlinear components
  epbasis = eigen(pbasis)  # eigen-decomposition
  upbasis = t(t(epbasis$vectors[, 1 : (nknots + 2)]) / sqrt(epbasis$values[1 : (nknots + 2)]))
  # concentration of linear and nonlinear parts to form new base
  sbasis = cbind(1, utime, crossprod(t(sbasis), upbasis))

  # calculate integral penalty
  pll = matrix(c(uknot - lknot, (uknot^2 - lknot^2) / 2, (uknot^2 - lknot^2)/2, (uknot^3 - lknot^3) / 3), nrow = 2)
  pln = tcrossprod(inprod(mobasis, bsbasis), t(upbasis)); pnn = tcrossprod(crossprod(upbasis, ibasis), t(upbasis))
  # combine all penalties into one matrix
  penat = cbind(rbind(pll, t(pln)), rbind(pln, pnn))

  list(lknot = lknot, uknot = uknot, iknot = iknot, sbasis = sbasis, penat = penat)
}

### Function to generate initial baseline coefficients

#' Title
#'
#' @param xobs
#' @param zobs
#' @param sbasis
#' @param svar
#' @param fscore
#' @param fpc
#' @param time
#'
#' @return
#' @export
#'
#' @examples
get_initial_zcoef = function(xobs, zobs, sbasis, svar, fscore, fpc, time) {
  utime = sort(unique(time)); numc = length(svar); numk = ncol(fscore) / numc; nums = ncol(sbasis); numz = ncol(zobs); numx = nrow(xobs)
  # get position of nodes
  cindex = c(1, which(time[-1] < time[-length(time)]) + 1, length(time) + 1); findex = seq(1, numc * numk + numc, by = numk)
  # initialize smoothness parameter
  eta = matrix(1e-8, numc, numz)

  zcoef = matrix(0, numc, numz * nums)
  # recursively sample for each node
  for(j in 1 : numc) {
    # match position for the current node
    cposition = cindex[j] : (cindex[j + 1] - 1); tposition = match(time[cposition], utime); fposition = findex[j] : (findex[j + 1] - 1)

    bsum = Bsum = 0
    for(i in 1 : numx) {
      nposition = which(!is.na(xobs[i, cposition])); ncposition = cposition[nposition]; ntposition = tposition[nposition]
      # remove the effect of principal components
      yobs = xobs[i, ncposition, drop = FALSE] - tcrossprod(fscore[i, fposition], fpc[[j]][nposition, , drop = FALSE])
      # sum the contribution of mean
      bsum = bsum + tcrossprod(yobs, kronecker(zobs[i, ], t(sbasis[ntposition, , drop = FALSE]))) / svar[j]
      # sum the contribution of precision
      Bsum = Bsum + tcrossprod(kronecker(zobs[i, ], t(sbasis[ntposition, , drop = FALSE]))) / svar[j]
    }

    # stack the penalty parameter
    pens = NULL; for(l in 1 : numz) pens = c(pens, rep(1e-8, 2), rep(eta[j, l], nums - 2))
    # Cholesky factorization of posterior precision matrix
    cholmat = chol(Bsum + diag(pens))

    # sample from the posterior distribution
    zcoef[j, ] = backsolve(cholmat, forwardsolve(t(cholmat), as.vector(bsum)) + rnorm(length(bsum)))

    # update the smoothness parameter
    for(l in 1 : numz) {
      shape = (nums - 2) / 2 + 0.001; rate = crossprod(zcoef[j, ((l - 1) * nums + 1) : (l * nums)][- c(1, 2)]) / 2 + 0.001
      # sample from the posterior gamma distribution
      eta[j, l] = rgamma(1, shape, rate)
    }
  }

  return(list(zcoef = zcoef, eta = eta))
}

### Function to generate initial principal components

get_initial_scoef = function(xobs, zobs, zcoef, sbasis, penat, svar, fscore, fpc, time) {
  utime = sort(unique(time)); numc = length(svar); numk = ncol(fscore) / numc; numx = nrow(xobs); nums = ncol(sbasis)
  # get positions of nodes
  cindex = c(1, which(time[-1] < time[-length(time)]) + 1, length(time) + 1); findex = seq(1, numc * numk + numc, by = numk)
  # initialize smoothness parameter
  lambda = rep(1e-8, numk)

  scoef = matrix(0, nums, numk)
  # recursively update for components
  for(k in 1 : numk) {
    bsum = Bsum = 0

    for(j in 1 : numc) {
      cposition = cindex[j] : (cindex[j + 1] - 1); tposition = match(time[cposition], utime); fposition = findex[j] : (findex[j + 1] - 1)
      # sum over non-missing values
      for(i in 1 : numx) {
        nposition = which(!is.na(xobs[i, cposition])); ncposition = cposition[nposition]; ntposition = tposition[nposition]
        # remove the baseline effects
        yobs = xobs[i, ncposition, drop = FALSE] - crossprod(zcoef[j, ], kronecker(zobs[i, ], t(sbasis[ntposition, , drop = FALSE])))
        # get contribution for mean
        bsum = bsum + crossprod(sbasis[ntposition, , drop = FALSE], t(yobs - tcrossprod(fscore[i, fposition[-k]], fpc[[j]][nposition, -k, drop = FALSE])) / svar[j] * fscore[i, fposition[k]])
        # get contribution from precision
        Bsum = Bsum + fscore[i, fposition[k]] * fscore[i, fposition[k]] / svar[j] * crossprod(sbasis[ntposition, , drop = FALSE])
      }
    }

    cholmat = chol(Bsum + diag(c(rep(1e-8, 2), rep(lambda[k], (nums - 2)))))
    # calculate the unconstrained vector
    ucont = backsolve(cholmat, forwardsolve(t(cholmat), bsum))
    # use sequential orthogonality for initialization
    if(k > 1){
      ocont = crossprod(t(penat), scoef[, 1 : (k - 1)])
      scont = backsolve(cholmat, forwardsolve(t(cholmat), ocont))
      ucont = ucont - crossprod(t(scont), crossprod(t(chol2inv(chol(crossprod(ocont, scont)))), crossprod(ocont, ucont)))
    }

    # get the final estimate by re-normalization
    unorm = sqrt(as.numeric(tcrossprod(crossprod(ucont, penat), t(ucont)))); scoef[, k] = ucont / unorm
    fscore[, seq(from = k, to = numc * numk, by = numk)] = fscore[, seq(from = k, to = numc * numk, by = numk)] * unorm

    # conditional maximum likelihood estimate of Lambda
    lambda[k] = (numc - 2) / crossprod(scoef[- c(1, 2), k])

    for(j in 1 : numc) {
      tposition = match(time[cindex[j] : (cindex[j + 1] - 1)], utime)
      fpc[[j]][, k] = crossprod(t(sbasis[tposition, ]), scoef[, k])
    }
  }

  list(scoef = scoef, lambda = lambda, fscore = fscore)
}

### Function to generate initial principal components
get_initial_scoef = function(xobs, zobs, zcoef, sbasis, penat, svar, fscore, fpc, time) {
  utime = sort(unique(time)); numc = length(svar); numk = ncol(fscore) / numc; numx = nrow(xobs); nums = ncol(sbasis)
  # get positions of nodes
  cindex = c(1, which(time[-1] < time[-length(time)]) + 1, length(time) + 1); findex = seq(1, numc * numk + numc, by = numk)
  # initialize smoothness parameter
  lambda = rep(1e-8, numk)

  scoef = matrix(0, nums, numk)
  # recursively update for components
  for(k in 1 : numk) {
    bsum = Bsum = 0

    for(j in 1 : numc) {
      cposition = cindex[j] : (cindex[j + 1] - 1); tposition = match(time[cposition], utime); fposition = findex[j] : (findex[j + 1] - 1)
      # sum over non-missing values
      for(i in 1 : numx) {
        nposition = which(!is.na(xobs[i, cposition])); ncposition = cposition[nposition]; ntposition = tposition[nposition]
        # remove the baseline effects
        yobs = xobs[i, ncposition, drop = FALSE] - crossprod(zcoef[j, ], kronecker(zobs[i, ], t(sbasis[ntposition, , drop = FALSE])))
        # get contribution for mean
        bsum = bsum + crossprod(sbasis[ntposition, , drop = FALSE], t(yobs - tcrossprod(fscore[i, fposition[-k]], fpc[[j]][nposition, -k, drop = FALSE])) / svar[j] * fscore[i, fposition[k]])
        # get contribution from precision
        Bsum = Bsum + fscore[i, fposition[k]] * fscore[i, fposition[k]] / svar[j] * crossprod(sbasis[ntposition, , drop = FALSE])
      }
    }

    cholmat = chol(Bsum + diag(c(rep(1e-8, 2), rep(lambda[k], (nums - 2)))))
    # calculate the unconstrained vector
    ucont = backsolve(cholmat, forwardsolve(t(cholmat), bsum))
    # use sequential orthogonality for initialization
    if(k > 1){
      ocont = crossprod(t(penat), scoef[, 1 : (k - 1)])
      scont = backsolve(cholmat, forwardsolve(t(cholmat), ocont))
      ucont = ucont - crossprod(t(scont), crossprod(t(chol2inv(chol(crossprod(ocont, scont)))), crossprod(ocont, ucont)))
    }

    # get the final estimate by re-normalization
    unorm = sqrt(as.numeric(tcrossprod(crossprod(ucont, penat), t(ucont)))); scoef[, k] = ucont / unorm
    fscore[, seq(from = k, to = numc * numk, by = numk)] = fscore[, seq(from = k, to = numc * numk, by = numk)] * unorm

    # conditional maximum likelihood estimate of Lambda
    lambda[k] = (numc - 2) / crossprod(scoef[- c(1, 2), k])

    for(j in 1 : numc) {
      tposition = match(time[cindex[j] : (cindex[j + 1] - 1)], utime)
      fpc[[j]][, k] = crossprod(t(sbasis[tposition, ]), scoef[, k])
    }
  }

  list(scoef = scoef, lambda = lambda, fscore = fscore)
}

