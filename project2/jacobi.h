#ifndef JACOBI_H
#define JACOBI_H

#include <iostream>
#include <armadillo>
#include <cmath>
#include <fstream>
#include <iomanip>
#include "time.h"

using namespace std;
using namespace arma;

void initialize(string& outfile_name, int& number_of_simulations,int& amplificationFactor, int& N, double& rhoMax,int& maxIterations, double& tolerance, string& armadillo, string& interactionRepulsion, double& omega, double& convergenceLimit, int argc, char** argv );
void jacobi(mat &A, colvec &eigenValues, double tolerance, int maxIterations, int N, int *counter, mat& v);
mat get_eigenvecs(mat a, mat v, colvec eigenValues, int N);
double findMaxNonDiagonalElement(mat &A, int *k, int *l, int N);
void rotate(mat &aMatrix, int k, int l, int N, mat& v);
void createEigenvalueVector( mat A, colvec &eigenValues, int N );
void createTridiagonalMatrix( mat &A, int N, double rhoMax, double rhoMin, double *h, string& interactionRepulsion, double omega);
void calculateError(colvec eigenValues, double *computedError);

#endif // JACOBI_H
