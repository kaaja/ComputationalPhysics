#ifndef SOLVER_H
#define SOLVER_H
#include "planet.h"

class Solver
{
private:
    Planet planet;
    int N;
    double step, time, x, y, vx, vy, r, pi, FourPi2, finalTime;

public:
    Solver();
    Solver(int N_, double T_);

    ~Solver() {} // Destructor
    void addPlanet(Planet planet_);

    void forwardEuler();
};

#endif // SOLVER_H
