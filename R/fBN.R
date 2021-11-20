### Loading required packages

library(fda)
library(splines2)
library(RGeode)

### Function to generate spline information

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
