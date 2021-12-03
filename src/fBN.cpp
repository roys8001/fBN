// [[Rcpp::depends(RcppArmadillo)]]
# include <RcppArmadillo.h>

using namespace arma;

// [[Rcpp::export]]
double rgammatr(double a, double b, arma::vec range, double n = 1) {
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
bool checkDAG(arma::mat G) {
  arma::mat cG = G; bool isDAG = TRUE;

  // find non-leaf nodes
  arma::uvec index = find(sum(cG, 0) != 0);
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
Rcpp::List get_zcoef(arma::mat xobs, arma::mat zobs, arma::mat scoef, arma::mat sbasis, arma::mat fscore, arma::mat eta, arma::vec svar, arma::uvec tindex, arma::uvec cindex, arma::uvec findex) {
  arma::uword numc = svar.n_elem, nums = sbasis.n_cols, numz = zobs.n_cols, numx = xobs.n_rows;
  arma::mat zcoef = zeros<mat>(numc, numz * nums);
  // recursively update for each row
  for(arma::uword j = 0; j < numc; j ++) {
    arma::uvec tposit = tindex(span(cindex(j), cindex(j + 1) - 1)); // position of available time points

    arma::vec msum = zeros<vec>(numz * nums); arma::mat vsum = zeros<mat>(numz * nums, numz * nums), feval = fscore.cols(span(findex(j), findex(j + 1) - 1)) * trans(sbasis.rows(tposit) * scoef);
    // sum over all samples (excluding missing values)
    for(arma::uword i = 0; i < numx; i ++) {
      arma::vec yobs = trans(xobs.row(i)); yobs = yobs(span(cindex(j), cindex(j + 1) - 1)) - trans(feval.row(i)); uvec uposit = find_finite(yobs);
      // calculate mean and covariance matrix
      msum = msum + trans(trans(yobs(uposit)) * kron(zobs.row(i), sbasis.rows(tposit(uposit))) / svar(j));
      vsum = vsum + trans(kron(zobs.row(i), sbasis.rows(tposit(uposit)))) * kron(zobs.row(i), sbasis.rows(tposit(uposit))) / svar(j);
    }

    arma::vec penat = zeros<vec>(numz * nums);
    // calculate spline penalization
    for(arma::uword k = 0; k < numz; k ++) penat(span(k * nums, (k + 1) * nums - 1)) = join_cols(1e-8 * ones<vec>(2), eta(j, k) * ones<vec>(nums - 2));
    // perform cholesky decomposition of covariance matrix
    arma::mat cholmat = chol(vsum + diagmat(penat));
    // generate from posterior normal distribution
    zcoef.row(j) = trans(solve(trimatu(cholmat), solve(trimatl(trans(cholmat)), msum) + randn(msum.n_elem)));

    // update penalty parameter
    arma::vec zsquare = trans(zcoef.row(j) % zcoef.row(j));
    for(arma::uword l = 0; l < numz; l ++) {
      double shape = (nums - 2) / 2 + 0.001, rate = accu(zsquare(span(l * nums + 2, (l + 1) * nums - 1))) / 2 + 0.001;
      // generate from posterior gamma distribution
      eta(j, l) = randg(distr_param(shape, 1 / rate));
    }
  }

  return Rcpp::List::create(Rcpp::Named("zcoef") = zcoef, Rcpp::Named("eta") = eta);
}

// [[Rcpp::export]]
Rcpp::List get_scoef(arma::mat xobs, arma::mat zobs, arma::mat zcoef, arma::mat scoef, arma::mat sbasis, arma::mat penat, arma::mat fscore, arma::vec svar, arma::vec lambda, arma::uvec tindex, arma::uvec cindex, arma::uvec findex, bool orderLambdas = 1) {
  arma::uword numc = svar.n_elem, numk = scoef.n_cols, nums = sbasis.n_cols, numx = xobs.n_rows;
  arma::uvec index = linspace<uvec>(0, numk - 1, numk);
  // recursively update for each column
  for(arma::uword k = 0; k < numk; k ++) {
    arma::vec msum = zeros<vec>(nums);
    arma::mat vsum = zeros<mat>(nums, nums);
    // sum over all nodes
    for(arma::uword j = 0; j < numc; j ++) {
      arma::uvec cposit = linspace<uvec>(cindex(j), cindex(j + 1) - 1, cindex(j + 1) - cindex(j));
      arma::uvec fposit = linspace<uvec>(findex(j), findex(j + 1) - 1, findex(j + 1) - findex(j));
      arma::uvec tposit = tindex(cposit); // position of available time points

      arma::mat feval = fscore.cols(fposit(find(index != k))) * trans(sbasis.rows(tposit) * scoef.cols(find(index != k)));
      // sum over all samples (excluding missing values)
      for(arma::uword i = 0; i < numx; i ++) {
        arma::vec yobs = trans(xobs.row(i)); yobs = yobs(cposit) - trans(zcoef.row(j) * trans(kron(zobs.row(i), sbasis.rows(tposit))));
        yobs = yobs - trans(feval.row(i)); uvec uposit = find_finite(yobs); vec ftemp = trans(fscore.row(i));
        // calculate mean and covariance matrix
        msum = msum + trans(ftemp(fposit(k)) * trans(yobs(uposit)) * sbasis.rows(tposit(uposit))) / svar(j);
        vsum = vsum + accu(ftemp(fposit(k)) * ftemp(fposit(k))) * trans(sbasis.rows(tposit(uposit))) * sbasis.rows(tposit(uposit)) / svar(j);
      }
    }

    // perform cholesky decomposition of covariance matrix
    arma::mat cholmat = arma::chol(vsum + diagmat(join_cols(1e-8 * ones<vec>(2), lambda(k) * ones<vec>(nums - 2))));
    // calculate joint orthogonality constraints
    arma::mat ocont = penat * scoef.cols(find(index != k));
    // sample the unconstrained vector
    arma::vec ucont = solve(trimatu(cholmat), solve(trimatl(trans(cholmat)), msum) + randn(msum.n_elem));
    // compute the vector of shift
    arma::mat scont = solve(trimatu(cholmat), solve(trimatl(trans(cholmat)), ocont));
    // calculate the constrained vector
    arma::vec vcont = ucont - scont * inv(trans(ocont) * scont) * trans(ocont) * ucont;

    scoef.col(k) = vcont / sqrt(accu(trans(vcont) * penat * vcont));

    // update penalty parameter
    arma::vec sctemp = scoef.col(k); double shape = (nums - 2) / 2 + 0.001, rate = accu(pow(sctemp(span(2, sctemp.n_elem - 1)), 2)) / 2 + 0.01;
    // set upper and lower bound constrain
    double lbound = 1e-8, ubound = 1e8; if(orderLambdas == 1) {
      if(k != 0) ubound = lambda(k - 1);
      if(k != (numk - 1)) lbound = lambda(k + 1);
    }

    arma::vec range = zeros<vec>(2); range(0) = lbound; range(1) = ubound;
    // generate lambda from truncated distribution
    lambda(k) = rgammatr(shape, rate, range);
  }

  return Rcpp::List::create(Rcpp::Named("scoef") = scoef, Rcpp::Named("lambda") = lambda);
}

// [[Rcpp::export]]
arma::mat get_fscore(arma::mat xobs, arma::mat zobs, arma::mat zcoef, arma::mat scoef, arma::mat sbasis, arma::mat fvar, arma::mat fcoef, arma::vec svar, arma::uvec tindex, arma::uvec cindex, arma::uvec findex) {
  arma::uword numx = xobs.n_rows, numc = svar.n_elem, numk = scoef.n_cols;
  arma::mat fscore = zeros<mat>(numx, numc * numk);

  // sum over all samples
  for(arma::uword i = 0; i < numx; i ++) {
    arma::mat var = zeros<mat>(numc * numk, numc * numk);
    arma::vec mean = zeros<vec>(numc * numk);
    // recursively for each node
    for(arma::uword j = 0; j < numc; j ++) {
      arma::uvec cposit = linspace<uvec>(cindex(j), cindex(j + 1) - 1, cindex(j + 1) - cindex(j));
      arma::uvec fposit = linspace<uvec>(findex(j), findex(j + 1) - 1, findex(j + 1) - findex(j));
      arma::uvec tposit = tindex(cposit); // position of available time points

      arma::vec yobs = trans(xobs.row(i));
      yobs = yobs(cposit) - trans(zcoef.row(j) * trans(kron(zobs.row(i), sbasis.rows(tposit))));
      arma::uvec uposit = find_finite(yobs);
      arma::mat ftemp = sbasis.rows(tposit(uposit)) * scoef;
      var(fposit, fposit) = trans(ftemp) * ftemp / svar(j);
      // calculate mean
      mean(fposit) = trans(trans(yobs(uposit)) * ftemp) / svar(j);
    }

    // calculate cholesky factorization
    arma::mat fcov = diagmat(1 / sqrt(fvar.row(i))) * (diagmat(ones<vec>(numc * numk)) - fcoef), cholmat = chol(var + trans(fcov) * fcov);
    // generate sample-specific score from posterior
    fscore.row(i) = trans(solve(trimatu(cholmat), solve(trimatl(trans(cholmat)), mean) + randn(mean.n_elem)));
  }

  return fscore;
}


// [[Rcpp::export]]
arma::vec get_svar(arma::mat xobs, arma::mat zobs, arma::mat zcoef, arma::mat scoef, arma::mat sbasis, arma::mat fscore, arma::uvec tindex, arma::uvec cindex, arma::uvec findex) {
  arma::uword numc = fscore.n_cols / scoef.n_cols, numx = xobs.n_rows;
  arma::vec svar = zeros<vec>(numc);
  // element-wise update
  for(arma::uword j = 0; j < numc; j ++) {
    arma::uvec cposit = linspace<uvec>(cindex(j), cindex(j + 1) - 1, cindex(j + 1) - cindex(j));
    arma::uvec fposit = linspace<uvec>(findex(j), findex(j + 1) - 1, findex(j + 1) - findex(j));
    arma::uvec tposit = tindex(cposit); // position of available time points
    // calculate shape parameter
    arma::mat xtemp = xobs.cols(cposit);
    xtemp = xtemp(find_finite(xtemp));
    double bnew = 0.01;
    double anew = 0.01 + xtemp.n_elem / 2;

    arma::mat feval = fscore.cols(fposit) * trans(sbasis.rows(tposit) * scoef);
    // sum over all samples (excluding missing values)
    for(arma::uword i = 0; i < numx; i ++) {
      arma::vec yobs = trans(xobs.row(i));
      yobs = yobs(cposit) - trans(zcoef.row(j) * trans(kron(zobs.row(i), sbasis.rows(tposit))));
      yobs = yobs - trans(feval.row(i)); uvec uposit = find_finite(yobs);
      // calculate rate parameter
      bnew = bnew + accu(yobs(uposit) % yobs(uposit)) / 2;
    }

    // sample from the posterior inverse gamma distribution
    svar(j) = 1 / randg(distr_param(anew, 1 / bnew));
  }

  return svar;
}


// [[Rcpp::export]]
Rcpp::List get_ginfo(arma::mat fscore, arma::mat fstructure, arma::mat gstructure, arma::mat fvar, arma::umat cmat, arma::vec degree, arma::vec mass, arma::uvec findex, double cvar, double eprob) {
  arma::uword numx = fscore.n_rows, numc = gstructure.n_cols, numk = fstructure.n_cols;

  // update graph structure
  for(arma::uword j = 0; j < numk; j ++) {
    for(arma::uword l = 0; l < numk; l ++) {
      arma::mat new_fstructure = fstructure, new_gsturecture = zeros<mat>(numc, numc);
      new_fstructure(j, l) = 1 - new_fstructure(j, l);
      // update outer structure
      for(arma::uword k = 0; k < numc; k ++) {
        for(arma::uword s = 0; s < numc; s ++) {
          arma::uvec row_fposit = linspace<uvec>(findex(k), findex(k + 1) - 1, findex(k + 1) - findex(k));
          arma::uvec col_fposit = linspace<uvec>(findex(s), findex(s + 1) - 1, findex(s + 1) - findex(s));

          if(accu(new_fstructure(row_fposit, col_fposit)) != 0) new_gsturecture(k, s) = 1;
        }
      }
      // check and update new structure
      if(checkDAG(new_gsturecture)) {
        arma::vec snum = fvar.col(j), sden = fvar.col(j);

        if(accu(fstructure.row(j)) != 0) sden = sden + cvar * sum(square(fscore.cols(find(fstructure.row(j)))), 1);
        if(accu(new_fstructure.row(j)) != 0) snum = snum + cvar * sum(square(fscore.cols(find(new_fstructure.row(j)))), 1);

        double den = accu(log_normpdf(fscore.col(j), zeros<vec>(numx), sqrt(sden))) + fstructure(j, l) * log(eprob) + (1 - fstructure(j, l)) * log(1 - eprob);
        double num = accu(log_normpdf(fscore.col(j), zeros<vec>(numx), sqrt(snum))) + new_fstructure(j, l) * log(eprob) + (1 - new_fstructure(j, l)) * log(1 - eprob);

        if((den - num) < log(1 / randu() - 1)) {fstructure = new_fstructure; gstructure = new_gsturecture;}
      }
    }
  }

  // update coefficient given graph structure
  arma::mat fcoef = zeros<mat>(numk, numk);
  // recursively for each coordinate
  for(arma::uword j = 0; j < numk; j ++) {
    arma::uvec oindex = find(fstructure.row(j) == 1);
    arma::uword ind = 0;

    if(oindex.n_elem != 0) {
      arma::mat ftemp = diagmat(1 / sqrt(fvar.col(j))) * fscore.cols(find(fstructure.row(j))), cholmat = chol(trans(ftemp) * ftemp + diagmat(ones<vec>(accu(fstructure.row(j)))) / cvar);
      arma::vec mean = trans(ftemp) * (fscore.col(j) / sqrt(fvar.col(j)));
      // calculate nonzero values
      arma::vec fvec = solve(trimatu(cholmat), solve(trimatl(trans(cholmat)), mean) + randn(mean.n_elem));

      for(arma::uword l = 0; l < numk; l ++) {
        if(any(oindex == l)) {fcoef(j, l) = fvec(ind); ind = ind + 1;} else fcoef(j, l) = 0;
      }
    } else fcoef.row(j) = trans(zeros<vec>(numk));
  }

  // update variance of nonzero coefficients
  double anew = 1 + accu(fstructure != 0) / 2, bnew = 1 + accu(pow(fcoef, 2)) / 2;
  // sample from the posterior
  cvar = 1 / randg(distr_param(anew, 1 / bnew));

  // update edge probability
  double cnew = 0.5 + accu(fstructure), dnew = 0.5 + accu(1 - fstructure) - numk, enew = 1;
  // sample from the posterior
  double apart = randg(distr_param(cnew, enew)), bpart = randg(distr_param(dnew, enew));
  eprob = apart / (apart + bpart);

  arma::mat cresidual = fscore - fscore * trans(fcoef);
  arma::uvec index = linspace<uvec>(0, numx - 1, numx);
  // update the Dirichlet process model
  for(arma::uword j = 0; j < numk; j ++) {
    for(arma::uword i = 0; i < numx; i ++) {
      arma::uvec ctemp = cmat.col(j);
      arma::uword ctotal = max(ctemp(find(index != i)));
      arma::vec cprob = zeros<vec>(ctotal + 1);
      // calculate probability for existing clusters
      for(arma::uword k = 1; k <= ctotal; k ++) {
        arma::vec btemp = cresidual.col(j);
        btemp = btemp(find(index != i));
        btemp = btemp(find(ctemp(find(index != i)) == k));
        double acount = accu(ctemp(find(index != i)) == k), bcount = accu(pow(btemp, 2));

        cprob(k - 1) = lgamma((degree(j) + acount + 1) / 2) - (degree(j) + acount + 1) / 2 * log((degree(j) + bcount + pow(cresidual(i, j), 2)) / 2);
        cprob(k - 1) = cprob(k - 1) + log(acount) - lgamma((degree(j) + acount) / 2) + (degree(j) + acount) / 2 * log((degree(j) + bcount) / 2);
      }
      // calculate probability of new cluster
      cprob(ctotal) = lgamma((degree(j) + 1) / 2) - (degree(j) + 1) / 2 * log((degree(j) + pow(cresidual(i, j), 2)) / 2);
      cprob(ctotal) = cprob(ctotal) + log(mass(j)) - lgamma(degree(j) / 2) + degree(j) / 2 * log(degree(j) / 2);
      // sample new cluster index
      cmat(i, j) = Rcpp::as<uword>(sample(Rcpp::Named("x") = linspace<uvec>(1, ctotal + 1, ctotal + 1), Rcpp::Named("size") = 1, Rcpp::Named("prob") = exp(cprob)));

      // remove empty clusters
      for(arma::uword k = 1; k <= (ctotal + 1); k ++) {
        if(accu(cmat.col(j) == k) == 0) {for(arma::uword t = 0; t < numx; t ++) {if(cmat(t, j) > k) cmat(t, j) = cmat(t, j) - 1;}}
      }
    }

    // update cluster specific parameter
    for(arma::uword k = 1; k <= max(cmat.col(j)); k ++) {
      arma::vec rtemp = cresidual.col(j);
      double anew = degree(j) + accu(cmat.col(j) == k) / 2, bnew = degree(j) + accu(pow(rtemp(find(cmat.col(j) == k)), 2)) / 2;
      // sample from the posterior
      double stemp = 1 / randg(distr_param(anew, 1 / bnew));
      for(arma::uword i = 0; i < numx; i ++) {if(cmat(i, j) == k) fvar(i, j) = stemp;}
    }

    // update mass parameter
    double cc = 1; double aphi = randg(distr_param(mass(j) + 1, cc)), bphi = randg(distr_param(numx, cc)), phi = aphi / (aphi + bphi);
    double pi = (1 + max(cmat.col(j)) - 1) / (numx * (1 - log(phi))); pi = pi / (1 + pi);
    // augmentation with pi
    if(randu() < pi) mass(j) = randg(distr_param(1 + max(cmat.col(j)), 1 / (1 - log(phi)))); else mass(j) = randg(distr_param(1 + max(cmat.col(j)) - 1, 1 / (1 - log(phi))));

    // update degree parameter: U(0, 10)
    double new_degree = 10 * randu();
    // calculate acceptance rate
    double den = accu(log(Rcpp::as<vec>(dgamma(Rcpp::Named("x") = 1 / unique(fvar.col(j)), Rcpp::Named("shape") = degree(j) / 2, Rcpp::Named("rate") = degree(j) / 2))));
    double num = accu(log(Rcpp::as<vec>(dgamma(Rcpp::Named("x") = 1 / unique(fvar.col(j)), Rcpp::Named("shape") = new_degree / 2, Rcpp::Named("rate") = new_degree / 2))));

    if((num - den) >= log(randu())) degree(j) = new_degree;
  }

  return Rcpp::List::create(Rcpp::Named("gstructure") = gstructure, Rcpp::Named("fstructure") = fstructure, Rcpp::Named("fcoef") = fcoef, Rcpp::Named("fvar") = fvar,
                                        Rcpp::Named("cmat") = cmat, Rcpp::Named("degree") = degree, Rcpp::Named("mass") = mass, Rcpp::Named("cvar") = cvar, Rcpp::Named("eprob") = eprob);
}

