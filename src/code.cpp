#include <RcppArmadillo.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

/*** R
# library(EloRating)
# library(microbenchmark)
*/

// [[Rcpp::export]]
arma::mat shuffle_mat(arma::mat inmat) {
  // randomize matrix
  arma::mat smat = inmat + trans(inmat) ;
  const int N = inmat.n_cols;
  arma::mat zmat = arma::zeros<arma::mat>(N, N);

  for (int i = 0; i < N; i++) {
    for (int j = (i + 1); j < N; j++) {
      arma::vec v = arma::randu<arma::vec>(1);
      arma::vec KK = round(smat(i, j) * v);
      zmat(i, j) = KK(0);
      zmat(j, i) = smat(j, i) - zmat(i, j) ;
    }
  }

  return zmat;
}

/*** R
# data("bonobos", package="EloRating")
# (x <- shuffle_mat(bonobos))
# sum(x)

*/


// [[Rcpp::export]]
arma::mat order_by_ds(arma::mat mat) {
  // first calc for original matrix
  // reorder matrix according to Pij DS
  arma::mat tmat = arma::trans(mat);
  arma::mat pmat = mat / (tmat + mat) ;
  int N = mat.n_cols;

  // replace NAs with 0s
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      if (std::isnan(pmat(j, i))) pmat(j, i) = 0;
    }
  }

  // get rowSums and colSums
  arma::colvec w = sum(pmat, 1);
  arma::rowvec l = sum(pmat, 0);
  arma::colvec w2 = pmat * w;
  arma::colvec l2 = trans(pmat) * trans(l);

  arma::rowvec DS = w.t() + w2.t() - l - l2.t();

  arma::uvec indices = sort_index(arma::conv_to<arma::vec>::from(DS), "descend");
  arma::mat outmat = mat(indices, indices);

  return outmat;
}


// [[Rcpp::export]]
arma::mat make_rank_matrix(arma::mat mat) {
  // create generic rank matrix
  int N2 = mat.n_cols;
  arma::mat DMAT(N2, N2);

  for(int i = 0; i < N2; i++) {
    DMAT(i, i) = 0;
    for(int j = (i + 1); j < N2; j++) {
      DMAT(i, j) = N2 * (i+1) - (N2*j-1) + N2*(j-1) + (j-i-1) - N2*i  ;
      DMAT(j ,i) = DMAT(i, j);
    }
  }
  return DMAT;
}

// [[Rcpp::export]]
double get_c(arma::mat inmat, arma::mat rankdistmat) {
  int N = inmat.n_cols;
  arma::mat PMAT = arma::zeros<arma::mat>(N, N);
  arma::mat DMAT = make_rank_matrix(inmat);
  for (int i = 0; i < N; i++) {
    for (int j = (i + 1); j < N; j++) {
      if( (inmat(i, j) + inmat(j, i)) > 0 ) {
        PMAT(i, j) = inmat(i, j) / (inmat(i, j) + inmat(j, i) * 1.0);
        PMAT(j, i) = 1.0 - PMAT(i, j);
      }
    }
  }

  arma::mat tPMAT = arma::trans(PMAT);
  arma::mat m1 = (PMAT - tPMAT) % DMAT;
  arma::mat m2 = (PMAT + tPMAT) % DMAT;
  double res = accu(trimatu(m1))/accu(trimatu(m2)) ;

  return res;
}



// [[Rcpp::export]]
List csteep(arma::mat inmat, int nrand) {
  // rank distance matrix
  arma::mat rmat = make_rank_matrix(inmat);
  // initialize objects
  NumericVector permres(nrand);
  double Cobs;
  double Hobs;
  double E;

  Cobs = get_c(order_by_ds(inmat), rmat);

  arma::mat testmat = inmat;
  for (int i = 0; i < nrand; i++) {
    testmat = order_by_ds(shuffle_mat(inmat));
    permres(i) = get_c(testmat, rmat);
  }


  E = mean(permres);

  if (Cobs >= E) {
    Hobs = (Cobs - E)/(1 - E);
  }
  if (Cobs < E) {
    Hobs = (Cobs - E) / E;
  }
  // List L = List::create(Named("name1") = v1 , _["name2"] = v2);
  return List::create(Named("Hobs") = Hobs, Named("Cobs") = Cobs, Named("E") = E, Named("randres") = permres);
}

/*** R
# data("bonobos", package="EloRating")
# csteep(bonobos, 10000)[1:3]
# csteepness::Hsteep(bonobos, 10000)[c(1, 4, 5)]
#
# csteep(devries98, 10000)[1:3]
# csteepness::Hsteep(devries98, 10000)[c(1, 4, 5)]
#
#
# EloRating::DS(bonobos, prop = "Pij")
*/
