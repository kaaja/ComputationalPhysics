#include "solver.h"
#include <fstream>
#include <iomanip>

ofstream ofile;

Solver:: Solver() { N = finalTime = 0; filename = "";}
Solver:: Solver(int N_, double finalTime_, string filename_)
{
    N = N_, finalTime = finalTime_;
    step = finalTime/double(N);
    time = 0.0;
    filename = filename_+string(".csv");
    ofile.open(filename);
    ofile << "time,x,y,vx/pi,vy/pi,potentialEnergy,kineticEnergy" << endl;

}

void Solver:: addPlanet(Planet planet_)
{
    planet = planet_;
}

void Solver:: forwardEuler()
{
    double potentialEnergy, kineticEnergy;
    x = planet.getInitialXPosition();
    y = planet.getInitialYPosition();
    vx = planet.getInitialXVelocity();
    vy = planet.getInitialYVelocity();
    r = sqrt(x*x + y*y);
    pi = acos(-1.0);
    FourPi2 = 4.*pi*pi;
    potentialEnergy = planet.getPotentialEnergy(r);
    kineticEnergy   = planet.getKineticEnergy(vx, vy);
    writeTofile(time, x, y, vx/pi, vy/pi, potentialEnergy, kineticEnergy);

    while (time < finalTime){
        x += step*vx;
        y += step*vy;
        vx -= step*FourPi2*x/(r*r*r);
        vy -= step*FourPi2*y/(r*r*r);
        r = sqrt(x*x + y*y);
        potentialEnergy = planet.getPotentialEnergy(r);
        kineticEnergy   = planet.getKineticEnergy(vx, vy);
        time += step;
        writeTofile(time, x, y, vx/pi, vy/pi, potentialEnergy, kineticEnergy);
    }
    ofile.close();
}

void Solver:: writeTofile(double time_, double x_, double y_, double vx_, double vy_, double potentialEnergy_, double kineticEnergy_ )
{
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    ofile << setw(15) << setprecision(8) << time_ << ", ";
    ofile << setw(15) << setprecision(8) << x_ << ", ";
    ofile << setw(15) << setprecision(8) << y_ << ", ";
    ofile << setw(15) << setprecision(8) << vx_ << ", ";
    ofile << setw(15) << setprecision(8) << vy_ << ", ";
    ofile << setw(15) << setprecision(8) << potentialEnergy_ << ", ";
    ofile << setw(15) << setprecision(8) << kineticEnergy_ << endl;
}

