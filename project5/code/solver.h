#ifndef SOLVER_H
#define SOLVER_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include "time.h"
#include <armadillo>

using namespace std;
using namespace arma;

class Solver
{
public:
    Solver( double dt_, double dx_, double theta_, double T_, int Nx_, int Nt_);
    void solve(string outfileName_);
private:
    string outfileName;
    double dt, dx, theta, T;
    int Nx, Nt;
    double alpha;
    double offDiagonalLhs;
    double diagonalLhs;
    double offDiagonalRhs;
    double diagonalRhs;
    //double * u, * computed_right_hand_side, *uOld;
    void generate_right_hand_side(double *computed_right_hand_side, double *uOld, double offDiagonalRhs, double diagonalRhs, int Nx);
    void gassianTridiagonalSymmetricSolver(double * computed_right_hand_side, double * u, double *uOld, double offDiagonalLhs, double diagonalLhs, int Nx);


};

#endif // SOLVER_H
