// 1st try on Jacobi algorithm for finding eigenvalues

#include <iostream>
#include <armadillo>
#include <cmath>


using namespace std;
using namespace arma;


void jacobi(mat A, colvec &eigenValues, double tolerance, int maxIterations, int N);
void findMaxNonDiagonalElement(mat &B, int *k, int *l, int N);
void computeB(mat aMatrix, mat &bMatrix, int k, int l);
void createEigenvalueVector( mat B, colvec &eigenValues, int N );



int main()
{
  int N = 3;
  double tolerance = pow(10, -5);
  int maxIterations = 100;

  mat A(3,3);
  A(0,0) = 1;
  A(0,1) = 2;
  A(0,2) = 3;
  A(1,0) = 2;
  A(1,1) = 1;
  A(1,2) = 3;
  A(2,0) = 5;
  A(2,1) = 1;
  A(2,2) = 3;


  colvec eigenValues  = zeros<colvec>(N);
  jacobi(A, eigenValues, tolerance, maxIterations, N);
  eigenValues.print("Eigenvalues =");
  return 0;
}

void jacobi(mat A, colvec &eigenValues, double tolerance, int maxIterations, int N){
    mat B;
    B = A;
    int k,l, counter;
    counter = 0;

    findMaxNonDiagonalElement(B, &k, &l, N);

    cout << "k: " << k << "l: " << l << endl;

    while ( B(k, l) > tolerance && counter < maxIterations ){
        findMaxNonDiagonalElement(B, &k, &l, N);
        computeB(A, B, k, l);
        counter += 1;
    }
    createEigenvalueVector( B, eigenValues,  N);
}

void computeB(mat aMatrix, mat &bMatrix, int k, int l){
    double tau, t, c, s;

    tau = (aMatrix(l,l) - aMatrix(k,k))/2*aMatrix(k,l);
    t = - tau + sqrt(1 + tau*tau);
    c = 1./sqrt(1+t*t);
    s = t*c;

    for (int i = 0; i < 2; i++){
        if ( (i != k) && (i != l) ) {
                bMatrix(i,i) = aMatrix(i, i);
                bMatrix(i,k) = aMatrix(i, k)*c - aMatrix(i, l)*s;
                bMatrix(i, l) = aMatrix(i, l)*c + aMatrix(i, k)*s;
        }
    bMatrix(k, k) = aMatrix(k, k)*c*c - 2.*aMatrix(k, l)*c*s + aMatrix(l, l)*s*s;
    bMatrix(l, l) = aMatrix(l, l)*c*c + 2.*aMatrix(k, l)*c*s + aMatrix(k, k)*s*s;
    bMatrix(k,l) = 0.0;
    bMatrix(l,k) = 0.0;
    }
}

void findMaxNonDiagonalElement(mat &B, int *k, int *l, int N){
    *k = 0;
    *l = 0;
    for (int row = 0; row < N; row++){
         for (int col = 1; col < N; col++){
             if( (row != col) && (B(row, col) > B(row, col -1)) ) {
                 *k = row;
                 *l = col;
             }
         }
    }
}

void createEigenvalueVector( mat B, colvec &eigenValues, int N ){
    for (int row = 0; row < N; row++){
        for (int col= 0; col < N; col++){
            if ( col == row )
                eigenValues(col) = B(row, col);
        }
    }
}

