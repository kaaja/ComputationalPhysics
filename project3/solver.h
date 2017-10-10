#ifndef SOLVER_H
#define SOLVER_H
#include "planet.h"
#include <string>
#include "time.h"
#include <fstream>
#include <iomanip>

class Solver
{
private:
    Planet planet;
    int N;
    double step, time, x, y, vx, vy, r, pi, FourPi2, finalTime, mass, potentialEnergy, kineticEnergy, angularMomentum, forceX, forceY, accelerationX, accelerationY, accelerationXOld, accelerationYOld, timeUsed;
    string filename;
    clock_t start, finish;

public:
    Solver();
    Solver(int N_, double T_, string filename_);

    ~Solver() {} // Destructor
    void addPlanet(Planet planet_);

    void forwardEuler();
    void velocityVerlet();
    void writeTofile(double time_, double x_, double y_, double vx_, double vy_, double potentialEnergy_, double kineticEnergy_, double angularMomentum_, double timeUsed_, double r_);


};

#endif // SOLVER_H
