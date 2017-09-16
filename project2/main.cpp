// 1st try on Jacobi algorithm for finding eigenvalues

#include <iostream>
#include <armadillo>
#include <cmath>


using namespace std;
using namespace arma;


void jacobi(mat &A, colvec &eigenValues, double tolerance, int maxIterations, int N);
double findMaxNonDiagonalElement(mat &A, int *k, int *l, int N);
void rotate(mat &aMatrix, int k, int l, int N);
void createEigenvalueVector( mat A, colvec &eigenValues, int N );



int main()
{
  int N = 3;
  double tolerance = pow(10, -10);
  int maxIterations = pow(10, 4);

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
  A.print("A:");
  return 0;
}

void jacobi(mat &A, colvec &eigenValues, double tolerance, int maxIterations, int N){
    int k,l, counter;
    counter = 0;

    double aMaxNonDiagonal = findMaxNonDiagonalElement(A, &k, &l, N);

    while ( fabs(aMaxNonDiagonal) > tolerance && counter < maxIterations ){
        aMaxNonDiagonal  = findMaxNonDiagonalElement(A, &k, &l, N);
        rotate(A, k, l, N);
        //findMaxNonDiagonalElement(A, &k, &l, N);
        counter += 1;
    }
    cout << "max a: " << fabs(A(k, l)) << endl;
    cout << "counter  " << counter << endl;
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
    //*k = 0;
    //*l = 0;
    double aMaxTemp = 0.0;
    for (int row = 0; row < N; row++){
         for (int col = row + 1; col < N; col++){
             if( row != col && fabs(A(row, col)) >  aMaxTemp ) {
                 *k = row;
                 *l = col;
                 aMaxTemp = fabs(A(row, col));
                 cout << "non diag max: " << aMaxTemp << endl;
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

