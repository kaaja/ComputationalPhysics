#ifndef EIGENVALUEBISECTION_H
#define EIGENVALUEBISECTION_H
#include "armadillo"
#include "float.h"

using namespace std;
using namespace arma;

colvec sturmSeq(mat &A, double &lam, int N);
colvec sturmSeqRevised(mat &A, double &lam, int N);
int numLambdas(colvec &p, int N);
int numLambdasRevised(colvec &p, int N);
colvec gerschgorin(mat &A, int N);
colvec lamRange(mat &A,int N, int oneForRevised);
double f(mat &A, double eigenvalueGuess, int N);
colvec eigenvals3(mat &A, int N, double bisectionAccuracy, int max_iterations, int oneForRevised);
double bisection(double (*func)(mat &A, double eigenvalueGuess, int N), double x1, double x2, double xacc, mat &A, int N, int max_iterations);
#endif // EIGENVALUEBISECTION_H
