#ifndef LANCZOS_H
#define LANCZOS_H

#include "armadillo"
#include <cmath>

using namespace std;
using namespace arma;

void lanczos(mat &A, colvec &alpha, colvec &beta, mat &Q, int N, int iterationNumber);
#endif // LANCZOS_H
