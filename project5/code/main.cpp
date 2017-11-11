#include <iostream>
#include "solver.h"

using namespace std;

int main()
{
    double dt = .1;
    double dx = .01;
    double theta = .5;
    double T = .1;
    Solver solution =  Solver(dt, dx, theta, T);
    solution.solve();
    return 0;
}
