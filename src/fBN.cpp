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

