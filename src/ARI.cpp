#include <Rcpp.h>
using namespace Rcpp;
//' @title calculate the value of adjusted Rand Index
//' @import Rcpp
//' @useDynLib SA24204166
//' @param k number of cluster
//' @param n nodes number
//' @param true_label the true label of nodes n x 1 vector
//' @param label the generated label of nodes n x 1 vector
//' @return the value of adjusted Rand Index
//' @export
//' @example
//' G<- DSBM(4,100,0.8,0.4)
//' H<-different_model(4,G$cg1,G$cg2,G$g1,G$D,100)
//' ARI(4,100,rep(1:4,each= 25),(H$km1$cluster))
// [[Rcpp::export]]
double ARI( int k,int n ,IntegerVector true_label, IntegerVector label) {
  NumericMatrix N(k, k);
  NumericVector a(k);
  NumericVector b(k);
  double index = 0;
  double aei = 0;
  double bei = 0;

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < k; ++j) {
      for (int l = 0; l < k; ++l) {
        if (true_label[i] - 1 == j && label[i] - 1 == l) {
          N(j, l)++;
        }
      }
    }
  }

  for (int i = 0; i < k; ++i) {
    a[i] = sum(N(i, _));
  }

  for (int i = 0; i < k; ++i) {
    b[i] = sum(N(_, i));
  }

  for (int i = 0; i < k; ++i) {
    for (int j = 0; j < k; ++j) {
      if (N(i, j) > 1) {
        index += N(i, j) * (N(i, j) - 1) / 2;
      }
    }
  }

  for (int i = 0; i < k; ++i) {
    if (a[i] > 1) {
      aei += a[i] * (a[i] - 1) / 2;
    }
  }

  for (int i = 0; i < k; ++i) {
    if (b[i] > 1) {
      bei += b[i] * (b[i] - 1) / 2;
    }
  }

  double ei = aei * bei * 2 / (n * (n - 1));
  double mri = (aei + bei) / 2;
  return (index - ei) / (mri - ei);
}
