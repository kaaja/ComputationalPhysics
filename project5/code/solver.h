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
    Solver(double dt_, double dx_, double theta_, double T_);
    void solve();
private:
    double dt, dx, theta, T;
    int Nx, Nt;
    double alpha;
    double offDiagonal;
    double diagonal;
    //double * u, * computed_right_hand_side, *uOld;
    void generate_right_hand_side(double *computed_right_hand_side, double *uOld, double offDiagonal, double diagonal, int Nx);
    void gassianTridiagonalSymmetricSolver(double * computed_right_hand_side, double * u, double *uOld, double offDiagonal, double diagonal, int Nx);


};

#endif // SOLVER_H
