#include <Rcpp.h>
using namespace Rcpp;

//' computes the no of elements that are in the ith group in a
//' and the jth group in b
//'
//' @param a,b integer vectors
//' @param i,j,n integers
//' @keywords internal
// [[Rcpp::export]]

double nij(NumericVector a, NumericVector b,int i,int j,int n) {
            int c;

            c=0;
            for (int l = 0; l < n; l++) {
            if(a(l) == i && b(l) == j){
            c++;
            }
            }
            return c;
}

//' computes Kendall's distance
//'
//' @param a,b integer vectors
//' @param n integer
//' @keywords internal
// [[Rcpp::export]]
  double kendR(NumericVector a, NumericVector b,int n) {
            double c;

            c=0;

            for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
            if(a(i)<a(j) && b(i) > b(j)){
            c++;
            }
            }
            }
            return c;
            }

