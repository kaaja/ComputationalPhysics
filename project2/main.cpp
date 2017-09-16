// 1st try on Jacobi algorithm for finding eigenvalues

#include <iostream>
#include <armadillo>
#include <cmath>

#define CATCH_CONFIG_RUNNER
#include "catch.hpp"

using namespace std;
using namespace arma;


void jacobi(mat &A, colvec &eigenValues, double tolerance, int maxIterations, int N);
double findMaxNonDiagonalElement(mat &A, int *k, int *l, int N);
void rotate(mat &aMatrix, int k, int l, int N);
void createEigenvalueVector( mat A, colvec &eigenValues, int N );
void createTridiagonalMatrix( mat &A, int N, double rhoMax, double rhoMin);


main(int argc, char* argv[]){ 
  int result = Catch::Session().run( argc, argv );
  int N = 3;
  double tolerance = pow(10, -10);
  int maxIterations = pow(10, 4);

  mat A;
  createTridiagonalMatrix(A, N, 3, 0);
  A.print("A:");


  colvec eigenValues  = zeros<colvec>(N);
  jacobi(A, eigenValues, tolerance, maxIterations, N);
  eigenValues.print("Eigenvalues =");
  return 0;
}

void jacobi(mat &A, colvec &eigenValues, double tolerance, int maxIterations, int N){
    int k,l, counter;
    counter = 0;

    double aMaxNonDiagonal = findMaxNonDiagonalElement(A, &k, &l, N);

    while ( fabs(aMaxNonDiagonal) > tolerance && counter < maxIterations ){
        aMaxNonDiagonal  = findMaxNonDiagonalElement(A, &k, &l, N);
        rotate(A, k, l, N);
        counter += 1;
    }
    createEigenvalueVector( A, eigenValues,  N);
}

void rotate(mat &aMatrix, int k, int l, int N){
    double c, s;
    if ( aMatrix(k,l) != 0.0 ) {
        double tau, t;
        tau = (aMatrix(l,l) - aMatrix(k,k))/(2.*aMatrix(k,l));
        if ( tau >= 0 ) {
          t = 1.0/(tau + sqrt(1.0 + tau*tau));
        }
        else {
          t = -1.0/(-tau +sqrt(1.0 + tau*tau));
        }

        c = 1/sqrt(1+t*t);
        s = c*t;
    }
    else {
        c = 1.0;
        s = 0.0;
    }

    double aIk, aIl, aKk, aLl;
    for (int i = 0; i < N; i++){
        if ( i != k && i != l ) {
                aIk = aMatrix(i, k);
                aIl = aMatrix(i, l);
                aMatrix(i,k) = aIk*c - aIl*s;
                aMatrix(k, i) = aMatrix(i, k);
                aMatrix(i, l) = aIl*c + aIk*s;
                aMatrix(l, i) = aMatrix(i, l);
        }
    }
    aKk = aMatrix(k,k);
    aLl = aMatrix(l, l);
    aMatrix(k, k) = aKk*c*c - 2.*aMatrix(k, l)*c*s + aLl*s*s ;
    aMatrix(l, l) = aLl*c*c + 2.*aMatrix(k, l)*c*s + aKk*s*s;
    aMatrix(k,l) = 0.0;
    aMatrix(l,k) = 0.0;
}

double findMaxNonDiagonalElement(mat &A, int *k, int *l, int N){
    double aMaxTemp = 0.0;
    for (int row = 0; row < N; row++){
         for (int col = row + 1; col < N; col++){
             if( row != col && fabs(A(row, col)) >  aMaxTemp ) {
                 *k = row;
                 *l = col;
                 aMaxTemp = fabs(A(row, col));
             }
         }
    }
    return aMaxTemp;
}

void createEigenvalueVector( mat A, colvec &eigenValues, int N ){
    for (int row = 0; row < N; row++){
        for (int col= 0; col < N; col++){
            if ( col == row )
                eigenValues(col) = A(row, col);
        }
    }
}

void createTridiagonalMatrix( mat &A, int N, double rhoMax, double rhoMin){
    double h = (rhoMax - rhoMin)/N;
    A.zeros(N,N);
    double offDiagonal = -1.0/(h*h);
    double diagonalFirstTerm = 2.0/(h*h);
    A(0, 1) =  offDiagonal;
    A(0, 0) = 2.0/(h*h);
    A(N-1, N-2) = offDiagonal;

    for (int row = 1; row < N-1; row++){
      A(row, row) = diagonalFirstTerm + pow(row*h,2);
      A(row, row - 1) = offDiagonal;
      A(row, row + 1) =  offDiagonal;
      }
    A(N-1, N-1) = 2.0/(h*h) + pow((N-1),2);
}

TEST_CASE( "Testing createTridiagonalMatrix", "[createTridiagonalMatrix]" ) {
  mat A;
  createTridiagonalMatrix(A, 10, 10, 0);
  REQUIRE( A(0,0)== 2);
}
