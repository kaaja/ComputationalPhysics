#ifndef LANCZOS2_H
#define LANCZOS2_H

#include "armadillo"
#include <cmath>
#include <string>

using namespace std;
using namespace arma;

void lanczos(colvec &eigenvalues, mat &A, colvec &alpha, colvec &beta, mat &Q, int N, int iterationNumber, string tridiag, string eigenvalueSolverm, int *stopIteration);
#endif // LANCZOS2_H
