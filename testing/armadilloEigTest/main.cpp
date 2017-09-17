
#include <iostream>
#include <armadillo>
#include <cmath>

using namespace std;
using namespace arma;

int main(){
  int dim = 5;
  mat A = randu<mat>(dim,dim);
  mat B = A.t()*A;  // generate a symmetric matrix

  vec eigval;

  eig_sym(eigval, B);
  eigval.print("Armadillo eigval: ");

}
