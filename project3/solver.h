#ifndef SOLVER_H
#define SOLVER_H
#include "planet.h"
#include "time.h"
#include <vector>

class Solver
{
private:
    int numberOfPlanets = 0;
    double step, x, y, vx, vy;
    double accelerationX = 0.0;
    double accelerationY =0.0;
    double accelerationXOld=0.0;
    double accelerationYOld = 0.0;
    vector<double> acc, accXVec, accYVec, accXVecOld, accYVecOld;
    string solverType;

    void forwardEuler(vector<Planet*> planets_, int numberOfplanets_, int iterationStart);
    void velocityVerlet(vector<Planet *> planets_, int numberOfplanets_, int iterationStart);


public:
    Solver(string solverType_);

    void solve(double step_,vector<Planet*> planets_, int numberOfplanets_, int iterationStart);

};

#endif // SOLVER_H
