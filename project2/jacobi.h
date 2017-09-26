#ifndef JACOBI_H
#define JACOBI_H

#include <iostream>
#include <armadillo>
#include <cmath>

using namespace std;
using namespace arma;


void jacobi(mat &A, colvec &eigenValues, double tolerance, int maxIterations, int N, int *counter, mat& v);
mat get_eigenvecs(mat a, mat v, colvec eigenValues, int N);
double findMaxNonDiagonalElement(mat &A, int *k, int *l, int N);
void rotate(mat &aMatrix, int k, int l, int N, mat& v);
void createEigenvalueVector( mat A, colvec &eigenValues, int N );

#endif // JACOBI_H
