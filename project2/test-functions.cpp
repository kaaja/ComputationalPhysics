#include "catch.hpp"
#include "jacobi.h"


TEST_CASE( "Testing maxDiagonal", "[findMaxNonDiagonalElement]" ) {
  // Checks that function finds max abs non-diagonal element in known symmetric matrix
  int k, l, N;
  N = 3;

  mat A(N,N);
  A(0,0) = 10.;
  A(0,1) = 1.;
  A(0,2) = -6.;
  A(1,0) = 1.;
  A(1,1) = 20.;
  A(1,2) = 3.;
  A(2,0) = -6.;
  A(2,1) = 3.;
  A(2,2) = 30.;

  double max = findMaxNonDiagonalElement(A, &k, &l, N);
  REQUIRE( max == 6.);
}


TEST_CASE( "Testing Jacobi eigenvalues", "[jacobi]" ) {
  // Checks that function finds max abs non-diagonal element in known symmetric matrix
  int k, l, N;
  N = 3;
  colvec eigenValues  = zeros<colvec>(N);
  vec EigvalArmadillo(N);
  double tolerance = pow(10,-10); // Tolerance for non-diagonals in Jacobi
  int maxIterations = pow(10,4);
  double tol = pow(10,-1); // Tolerance for difference jacobi and Armadillo
  mat v = eye<mat>(N,N);

  mat A(N,N);
  A(0,0) = 10.;
  A(0,1) = 1.;
  A(0,2) = -6.;
  A(1,0) = 1.;
  A(1,1) = 20.;
  A(1,2) = 3.;
  A(2,0) = -6.;
  A(2,1) = 3.;
  A(2,2) = 30.;
  int counter;
  eig_sym(EigvalArmadillo, A);
  jacobi(A, eigenValues, tolerance, maxIterations, N, &counter, v);
  for (int i = 0; i < N; i++){
          REQUIRE( fabs((eigenValues(i) - EigvalArmadillo(i))/EigvalArmadillo(i)) < tol);
  }
}




/*
TEST_CASE( "Testing jacobi", "[jacobi]" ) {
  mat A;
  int N = pow(10,3)
  createTridiagonalMatrix(A, N, 10, 0);
  colvec eigenValues;
  jacobi(A, eigenValues, pow(10,-7), pow(10,3), N );
  vec Eigval(N);
  eig_sym(Eigval, A);
  REQUIRE( A(0,0)== 2);
}
*/
