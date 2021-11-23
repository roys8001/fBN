// [[Rcpp::depends(RcppArmadillo)]]
# include <RcppArmadillo.h>

using namespace arma;

// [[Rcpp::export]]
double rgammatr(double a, double b, vec range, double n = 1) {
  // obtain environment containing function
  Rcpp::Environment RGeode = Rcpp::Environment::namespace_env("RGeode");
  // make function callable from Cpp
  Rcpp::Function trgamma = RGeode["rgammatr"];

  // return object
  return Rcpp::as<double>(trgamma(n, a, b, range));
}

Rcpp::Function dgamma("dgamma");
Rcpp::Function sample("sample");

// [[Rcpp::export]]
bool checkDAG(mat G) {
  mat cG = G; bool isDAG = TRUE;

  // find non-leaf nodes
  uvec index = find(sum(cG, 0) != 0);
  while(index.n_elem < cG.n_rows) {
    // delete leaf nodes
    cG = cG(index, index);
    if(accu(cG) == 0) break; else {
      // find remaining non-leaf nodes
      index = find(sum(cG, 0) != 0);
    }
  }

  if(accu(cG) != 0) isDAG = FALSE;

  return(isDAG);
}

// [[Rcpp::export]]
Rcpp::List get_zcoef(mat xobs, mat zobs, mat scoef, mat sbasis, mat fscore, mat eta, vec svar, uvec tindex, uvec cindex, uvec findex) {
  uword numc = svar.n_elem, nums = sbasis.n_cols, numz = zobs.n_cols, numx = xobs.n_rows; mat zcoef = zeros<mat>(numc, numz * nums);
  // recursively update for each row
  for(uword j = 0; j < numc; j ++) {
    uvec tposit = tindex(span(cindex(j), cindex(j + 1) - 1)); // position of available time points

    vec msum = zeros<vec>(numz * nums); mat vsum = zeros<mat>(numz * nums, numz * nums), feval = fscore.cols(span(findex(j), findex(j + 1) - 1)) * trans(sbasis.rows(tposit) * scoef);
    // sum over all samples (excluding missing values)
    for(uword i = 0; i < numx; i ++) {
      vec yobs = trans(xobs.row(i)); yobs = yobs(span(cindex(j), cindex(j + 1) - 1)) - trans(feval.row(i)); uvec uposit = find_finite(yobs);
      // calculate mean and covariance matrix
      msum = msum + trans(trans(yobs(uposit)) * kron(zobs.row(i), sbasis.rows(tposit(uposit))) / svar(j));
      vsum = vsum + trans(kron(zobs.row(i), sbasis.rows(tposit(uposit)))) * kron(zobs.row(i), sbasis.rows(tposit(uposit))) / svar(j);
    }

    vec penat = zeros<vec>(numz * nums);
    // calculate spline penalization
    for(uword k = 0; k < numz; k ++) penat(span(k * nums, (k + 1) * nums - 1)) = join_cols(1e-8 * ones<vec>(2), eta(j, k) * ones<vec>(nums - 2));
    // perform cholesky decomposition of covariance matrix
    mat cholmat = chol(vsum + diagmat(penat));
    // generate from posterior normal distribution
    zcoef.row(j) = trans(solve(trimatu(cholmat), solve(trimatl(trans(cholmat)), msum) + randn(msum.n_elem)));

    // update penalty parameter
    vec zsquare = trans(zcoef.row(j) % zcoef.row(j));
    for(uword l = 0; l < numz; l ++) {
      double shape = (nums - 2) / 2 + 0.001, rate = accu(zsquare(span(l * nums + 2, (l + 1) * nums - 1))) / 2 + 0.001;
      // generate from posterior gamma distribution
      eta(j, l) = randg(distr_param(shape, 1 / rate));
    }
  }

  return Rcpp::List::create(Rcpp::Named("zcoef") = zcoef, Rcpp::Named("eta") = eta);
}

// [[Rcpp::export]]
Rcpp::List get_scoef(mat xobs, mat zobs, mat zcoef, mat scoef, mat sbasis, mat penat, mat fscore, vec svar, vec lambda, uvec tindex, uvec cindex, uvec findex, bool orderLambdas = 1) {
  uword numc = svar.n_elem, numk = scoef.n_cols, nums = sbasis.n_cols, numx = xobs.n_rows; uvec index = linspace<uvec>(0, numk - 1, numk);
  // recursively update for each column
  for(uword k = 0; k < numk; k ++) {
    vec msum = zeros<vec>(nums); mat vsum = zeros<mat>(nums, nums);
    // sum over all nodes
    for(uword j = 0; j < numc; j ++) {
      uvec cposit = linspace<uvec>(cindex(j), cindex(j + 1) - 1, cindex(j + 1) - cindex(j));
      uvec fposit = linspace<uvec>(findex(j), findex(j + 1) - 1, findex(j + 1) - findex(j));
      uvec tposit = tindex(cposit); // position of available time points

      mat feval = fscore.cols(fposit(find(index != k))) * trans(sbasis.rows(tposit) * scoef.cols(find(index != k)));
      // sum over all samples (excluding missing values)
      for(uword i = 0; i < numx; i ++) {
        vec yobs = trans(xobs.row(i)); yobs = yobs(cposit) - trans(zcoef.row(j) * trans(kron(zobs.row(i), sbasis.rows(tposit))));
        yobs = yobs - trans(feval.row(i)); uvec uposit = find_finite(yobs); vec ftemp = trans(fscore.row(i));
        // calculate mean and covariance matrix
        msum = msum + trans(ftemp(fposit(k)) * trans(yobs(uposit)) * sbasis.rows(tposit(uposit))) / svar(j);
        vsum = vsum + accu(ftemp(fposit(k)) * ftemp(fposit(k))) * trans(sbasis.rows(tposit(uposit))) * sbasis.rows(tposit(uposit)) / svar(j);
      }
    }

    // perform cholesky decomposition of covariance matrix
    mat cholmat = chol(vsum + diagmat(join_cols(1e-8 * ones<vec>(2), lambda(k) * ones<vec>(nums - 2))));
    // calculate joint orthogonality constraints
    mat ocont = penat * scoef.cols(find(index != k));
    // sample the unconstrained vector
    vec ucont = solve(trimatu(cholmat), solve(trimatl(trans(cholmat)), msum) + randn(msum.n_elem));
    // compute the vector of shift
    mat scont = solve(trimatu(cholmat), solve(trimatl(trans(cholmat)), ocont));
    // calculate the constrained vector
    vec vcont = ucont - scont * inv(trans(ocont) * scont) * trans(ocont) * ucont;

    scoef.col(k) = vcont / sqrt(accu(trans(vcont) * penat * vcont));

    // update penalty parameter
    vec sctemp = scoef.col(k); double shape = (nums - 2) / 2 + 0.001, rate = accu(pow(sctemp(span(2, sctemp.n_elem - 1)), 2)) / 2 + 0.01;
    // set upper and lower bound constrain
    double lbound = 1e-8, ubound = 1e8; if(orderLambdas == 1) {
      if(k != 0) ubound = lambda(k - 1);
      if(k != (numk - 1)) lbound = lambda(k + 1);
    }

    vec range = zeros<vec>(2); range(0) = lbound; range(1) = ubound;
    // generate lambda from truncated distribution
    lambda(k) = rgammatr(shape, rate, range);
  }

  return Rcpp::List::create(Rcpp::Named("scoef") = scoef, Rcpp::Named("lambda") = lambda);
}

