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
    vector<double> accelerationsX;
    vector<double> accelerationsY;
    vector<double> acc, accXVec, accYVec, accXVecOld, accYVecOld;
    int N;
    int numberOfPlanets = 0;
    double step, time, x, y, vx, vy, r, finalTime, timeUsed;
    double accelerationX = 0.0;
    double accelerationY =0.0;
    double accelerationXOld=0.0;
    double accelerationYOld = 0.0;
    string filename;
    bool centerOfMassSystem;
    clock_t start, finish;

    double getCenterOfMassX();
    double getCenterOfMassY();
    void setSunVelocity();


public:
    Solver();
    Solver(int N_, double T_);

    ~Solver() {} // Destructor
    void addPlanet(Planet &planet_);

    void forwardEuler();
    void forwardEulerOld();
    void velocityVerlet();
    void velocityVerletOld();
    void alternativeForceVelocityVerlet(double beta_);

    void changeToCenterOfMassSystem();

};

#endif // SOLVER_H
