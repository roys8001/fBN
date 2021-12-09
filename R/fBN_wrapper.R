#' Title Functional Bayesian Network
#'
#' @param xobs a \enq{n × |tobs|} matrix where n is the number of samples. Its (i, j)the element indicates the discrete observation for functional covariate k at time point \enq{t^{h}_k} of sample i if \enq{j = \sum_{l = 1}^{k - 1}|t_l| + h} according to the arrangement of tobs. Missing values are allowed and specified as NaN.
#' @param zobs a \enq{n × q} matrix where q is the number of baseline scalar covariates
#' @param time zzzz
#' @param numk zzxxxx
#' @param nknots jjhg
#' @param tolCPV kdodjk
#' @param useAT higag
#'
#' @return
#' @export
#'
#' @examples
#' library(rmutil); library(orthogonalsplinebasis); library(igraph); library(splines)

#' numx = 500 # Number of sample points
#' numc = 10 # Number of variables
#' numz = 2 # Number of baseline scalar covariates
#' numk = 4 # Number of principal components
#' numt = 200 # Total number of time points
#' numxt = 10 # Number of timepoints for (sample, variable)
#' vstrength = 1
#' vnoise = 0.5

#' # simulation truth of outer graph structure: DAG
#' G0 = matrix(0, numc, numc);
#' G0[lower.tri(G0)] = as.numeric(runif(sum(lower.tri(G0))) < 2 / numc)
#' # simulation truth of inner graph structure
#' GS = matrix(0, numc * numk, numc * numk)
#' for(j in 2 : numc) {
#'   for(l in 1 : (j - 1)) {
#'     if(G0[j, l] == 1) {
#'       cindex = ((j - 1) * numk + 1) : (j * numk);
#'       findex = ((l - 1) * numk + 1) : (l * numk)
#'       # selection of edges
#'       GS[cindex, findex] = as.numeric(runif(numk * numk) < 0.5)
#'       while(sum(GS[cindex, findex]) == 0) GS[cindex, findex] = as.numeric(runif(numk * #' numk) < 0.5)
#'     }
#'   }
#' }

#' # nonzero signal strength
#' BS = matrix(0, numc * numk, numc * numk)
#' BS[GS != 0] = runif(sum(GS != 0), vstrength / 2, vstrength) # uniform
#' BS[GS != 0] = BS[GS != 0] * (1 - 2 * (runif(sum(GS != 0)) < 0.5))
#'
#' # generate principal scores
#' score = matrix(rlaplace(numx * numc * numk, s = vnoise/2), numx, numc * numk)
#' score = t(tcrossprod(solve(diag(numc * numk) - BS), score))
#'
#' # generate all available observation time
#' t0 = sort(runif(numt))
#' tx = array(0, dim = c(numx, numc, numxt))
#' # generate observation time for each sample and node
#' for(i in 1 : numx) for(j in 1 : numc) tx[i, j, ] = sample(t0, numxt)
#' # summarize observation time for each node
#' tobs = NULL
#' for(j in 1 : numc) tobs = c(tobs, sort(unique(as.vector(tx[, j, ]))))
#'
#' # generate functional principal components
#' lb = min(t0)
#' ub = max(t0)
#' ib = quantile(t0, c(0, 1))
#' fPC = evaluate(OBasis(knots = c(rep(lb, numk - 1), ib, rep(ub, numk - 1)), order = numk), t0)
#'
#' # generate baseline variables
#' zobs = matrix(rnorm(numx * numz), numx, numz)
#'
#' # polynomial basis for baseline effects
#' zbase = bs(t0);
#' nums = ncol(zbase);
#' ecoef = matrix(runif(numc * numz * nums, vstrength / 2, vstrength), numc, numz * nums)
#' ecoef = ecoef * (1 - 2 * (runif(length(ecoef)) < 0.5))
#'
#' # generate observation
#' xobs = matrix(0, numx, length(tobs))
#' cindex = c(1, which(tobs[-1] < tobs[-length(tobs)]) + 1, length(tobs) + 1)
#' e0 = rep(0, numc)
#' for(j in 1 : numc) {
#'   # node index
#'   fposition = ((j - 1) * numk + 1) : (j * numk); cposition = cindex[j] : (cindex[j + 1] - 1)
#'   # calculate signal strength
#'   ee = tcrossprod(score[, fposition], fPC[match(tobs[cposition], t0), ])
#'   for(l in 1 : numz) ee = ee + crossprod(t(zobs[, l]), tcrossprod(ecoef[j, ((l - 1) * nums + 1) : (l * nums)], zbase[match(tobs[cposition], t0), ]))
#'
#'   e0[j] = vnoise * mean(abs(ee))
#'   # match sample with observation time
#'   for(i in 1 : numx) {
#'     xposition = match(tx[i, j, ], tobs[cposition]); tposition = match(tx[i, j, ], t0)
#'     # generate observed data
#'     xobs[i, cposition[xposition]] = e0[j] * rnorm(numxt) + crossprod(ecoef[j, ], kronecker(zobs[i, ], t(zbase[tposition, ]))) + tcrossprod(score[i, fposition], fPC[tposition, ])
#'   }
#' }
#' # fill zeros with NA
#' xobs[xobs == 0] = NaN
#' get_initial_parameter(xobs = xobs, zobs = zobs, time = tobs)

get_initial_parameter = function(xobs, zobs, time, numk = NULL, nknots = 20, tolCPV = 0.9, useAT = FALSE) {
  utime = sort(unique(time)); cindex = c(1, which(time[-1] < time[-length(time)]) + 1, length(time) + 1)
  numc = length(cindex) - 1; numx = nrow(xobs)
  # get relevant spline information
  sinfo = get_sinfo(time, nknots); sbasis = sinfo$sbasis; penat = sinfo$penat

  # specify the time points to interpolate for initialization using SVD
  if(useAT) t = utime else t = c(sinfo$lknot, sinfo$iknot, sinfo$uknot)

  # compute a complete matrix without missing values (combined all samples and nodes)
  xc = array(0, c(numc * numx, length(t)))
  # iterate over each node
  for(j in 1 : numc) {
    yobs = xobs[, cindex[j] : (cindex[j + 1] - 1)]; scase = which(rowSums(!is.na(yobs)) != 0)
    # smooth across time if observation is partially missing
    yobs[scase, ] = t(apply(yobs[scase, ], 1, function(x) splinefun(time[cindex[j] : (cindex[j + 1] - 1)], x, method = 'natural') (time[cindex[j] : (cindex[j + 1] - 1)])))
    # smooth across sample if observation is completely missing
    yobs = apply(yobs, 2, function(x) {splinefun(1 : numx, x, method = 'natural') (1 : numx)})

    # xc is a stacked data matrix with numc time numx rows and |t| columns
    yobs = t(apply(yobs, 1, function(x) splinefun(time[cindex[j] : (cindex[j + 1] - 1)], x, method = 'natural') (t)))
    # remove the baseline effect
    xc[(1 : numx) + (j * numx - numx), ] = apply(yobs, 2, function(x) lm(x ~ zobs)$residual)
  }

  # compute SVD of the imputed data matrix and cumulative sums of the squares proportions
  svdx = svd(xc); cumv = cumsum(svdx$d^2 / sum(svdx$d^2))
  # select number of principal components based on cumulative sums if numk is unspecified
  if(is.null(numk)) numk = max(2, which(cumv >= tolCPV)[1])

  print(t(matrix(cumv, dimnames = list(paste('numk =', 1 : length(svdx$d))))))

  # get initial principal components and scores
  findex = seq(1, numc * numk + numc, by = numk); fpc = vector("list", numc); fscore = matrix(0, numx, numc * numk)
  # iterate over all nodes
  for(j in 1 : numc) {
    # if use all time points then subset for each node; otherwise smooth the SVD estimates and interpolate
    if(useAT) fpc[[j]] = as.matrix(svdx$v[match(time[cindex[j] : (cindex[j + 1] - 1)], utime), 1 : numk])
    if(!useAT) fpc[[j]] = apply(as.matrix(svdx$v[, 1 : numk]), 2, function(x) {splinefun(t, x, method = 'natural') (time[cindex[j] : (cindex[j + 1] - 1)])})
    # calculate principal scores
    fscore[, findex[j] : (findex[j + 1] - 1)] = (crossprod(t(svdx$u), diag(svdx$d)))[(1 : numx) + (j - 1) * numx, 1 : numk]
  }

  # initial sampling variance
  svar = rep(1, numc)
  # initial baseline coefficients
  zinfo = get_initial_zcoef(xobs, zobs, sbasis, svar, fscore, fpc, time)
  zcoef = zinfo$zcoef; eta = zinfo$eta
  # initial basis coefficients
  sinfo = get_initial_scoef(xobs, zobs, zcoef, sbasis, penat, svar, fscore, fpc, time)
  scoef = sinfo$scoef; lambda = sinfo$lambda; fscore = sinfo$fscore

  # initialize graph information (empty)
  gstructure = matrix(0, numc, numc); fstructure = matrix(0, numc * numk, numc * numk); mass = rep(1, numc * numk)
  fcoef = matrix(0, numc * numk, numc * numk); fvar = matrix(1, numx, numc * numk); cvar = 1
  cmat = matrix(1, numx, numc * numk); degree = rep(1, numc * numk); eprob = 0.5

  tindex = NULL; for(j in 1 : numc) tindex = c(tindex, match(time[cindex[j] : (cindex[j + 1] - 1)], sort(unique(time))))
  # match subscript for c++ input
  cindex = cindex - 1; findex = findex - 1; tindex = tindex - 1

  # run 10 iterations to obtain ordering
  for(iter in 1 : 10) {
    svar = get_svar(xobs, zobs, zcoef, scoef, sbasis, fscore, tindex, cindex, findex)
    # baseline effects
    zinfo = get_zcoef(xobs, zobs, scoef, sbasis, fscore, eta, svar, tindex, cindex, findex)
    zcoef = zinfo$zcoef; eta = zinfo$eta
    # spline information
    sinfo = get_scoef(xobs, zobs, zcoef, scoef, sbasis, penat, fscore, svar, lambda, tindex, cindex, findex, orderLambdas = 0)
    scoef = sinfo$scoef; lambda = sinfo$lambda
    # principal scores
    fscore = get_fscore(xobs, zobs, zcoef, scoef, sbasis, fvar, fcoef, svar, tindex, cindex, findex)
    # causal structure information
    ginfo = get_ginfo(fscore, fstructure, gstructure, fvar, cmat, degree, mass, findex, cvar, eprob)

    fvar = ginfo$fvar; fcoef = ginfo$fcoef; gstructure = ginfo$gstructure; fstructure = ginfo$fstructure
    cvar = ginfo$cvar; eprob = ginfo$eprob; cmat = ginfo$cmat; degree = ginfo$degree; mass = ginfo$mass
  }

  # order the k components based on the smoothing parameters lambda
  sorder = order(lambda, decreasing = TRUE); lambda = lambda[sorder]; scoef = scoef[, sorder]
  for(j in 1 : numc) {fposition = (findex[j] + 1) : findex[j + 1]; fscore[, fposition] = fscore[, fposition][, sorder]}

  # make sure the l = 1 curve is positive and might as well initialize all curves to have positive sums
  for(k in 1 : numk) {
    if(sum(crossprod(t(sbasis), scoef[, k])) < 0) {
      scoef[, k] = - scoef[, k]; fscore[, seq(from = k, to = numc * numk , by = numk)] = - fscore[, seq(from = k, to = numc * numk, by = numk)]
    }
  }

  list(fscore = fscore, zcoef = zcoef, eta = eta, scoef = scoef, lambda = lambda, svar = svar, fvar = fvar, fcoef = fcoef, numk = numk, penat = penat, findex = findex,
       gstructure = gstructure, fstructure = fstructure, cvar = cvar, eprob = eprob, cmat = cmat, degree = degree, mass = mass, sbasis = sbasis, cindex = cindex, tindex = tindex)
}
