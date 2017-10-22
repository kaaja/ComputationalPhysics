#ifndef System_H
#define System_H
#include "planet.h"
#include "solver.h"
#include "time.h"
#include <vector>

class System
{
private:
    Planet* planet;
    vector<Planet*> planets;
    int N;
    int numberOfPlanets = 0;
    double step, time, x, y, vx, vy, r, finalTime, timeUsed;

    string filename;
    bool centerOfMassSystem;
    clock_t start, finish;

    double getCenterOfMassX();
    double getCenterOfMassY();
    void setSunVelocity();


public:
    System();

    ~System() {} // Destructor
    void addPlanet(Planet &planet_);

    void simulate(Solver &solver, int N_, double finalTime_);

    void changeToCenterOfMassSystem();

};

#endif // System_H
