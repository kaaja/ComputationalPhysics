#ifndef SOLVER_H
#define SOLVER_H
#include "planet.h"
#include "time.h"
//#include <fstream>
//#include <iomanip>
#include <vector>

class Solver
{
private:
    Planet* planet;
    vector<Planet*> planets;
    int N;
    int numberOfPlanets = 0;
    double step, time, x, y, vx, vy, r, pi, FourPi2, finalTime, mass, potentialEnergy, kineticEnergy, angularMomentum, forceX, forceY, accelerationX, accelerationY, accelerationXOld, accelerationYOld, timeUsed;
    string filename, centerOfMassSystem;
    clock_t start, finish;

    double getCenterOfMassX();
    double getCenterOfMassY();
    void setSunVelocity();


public:
    Solver();
    Solver(int N_, double T_, string filename_);

    ~Solver() {} // Destructor
    void addPlanet(Planet &planet_);

    void forwardEuler();
    void velocityVerlet();
    void alternativeForceVelocityVerlet(double beta_);

    void changeToCenterOfMassSystem();

};

#endif // SOLVER_H
