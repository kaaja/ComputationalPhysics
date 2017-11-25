#ifndef TWODIMENSIONALDIFFUSIONSOLVER_H
#define TWODIMENSIONALDIFFUSIONSOLVER_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include "time.h"
#include <armadillo>
#include <omp.h>


using namespace std;
using namespace arma;

class TwoDimensionalDiffusionSolver
{
private:
    string outfileName;
    double dt, dx, dy, theta, T, alpha;
    int Nx, Ny, Nt;

    // Steady state solution
    void explicitScheme( double ** u, double ** uOld, int Nx, int Ny);
    void backwardEulerJacobi(double **u, double **uOld, int Nx, int Ny, int maxIterations, double maxDifference, double &iterationNumber);
    void backwardEulerGaussSeidel(double **u, double **uOld, int Nx, int Ny, int maxIterations, double maxDifference, double &iterationNumber);
    double ** CreateMatrix(int m, int n);
    void DestroyMatrix(double ** mat, int m, int n);
    void setMatrixAEqualMatrixB(double ** matrixA, double ** matrixB, int m, int n);

public:
    TwoDimensionalDiffusionSolver(double dt_, double dx_, double dy_, double theta_, double T_, int Nx_, int Ny_, int Nt_);
    void solve(string outfileName_, string method_, int thread_num_);
    double uSteadyState(double x, double y);

};

#endif // TWODIMENSIONALDIFFUSIONSOLVER_H
