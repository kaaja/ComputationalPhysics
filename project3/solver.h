#ifndef SOLVER_H
#define SOLVER_H
#include "planet.h"
#include <string>

class Solver
{
private:
    Planet planet;
    int N;
    double step, time, x, y, vx, vy, r, pi, FourPi2, finalTime;
    string filename;

public:
    Solver();
    Solver(int N_, double T_, string filename_);

    ~Solver() {} // Destructor
    void addPlanet(Planet planet_);

    void forwardEuler();
    void writeTofile(double time_, double x_, double y_, double vx_, double vy_ );
};

#endif // SOLVER_H
