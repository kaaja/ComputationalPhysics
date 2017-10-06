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
    filename = filename_;
    ofile.open(filename);
    ofile << ",time,x,y,vx/pi,vy/pi" << endl;

}

void Solver:: addPlanet(Planet planet_)
{
    planet = planet_;
}

void Solver:: forwardEuler()
{
    x = planet.getInitialXPosition();
    y = planet.getInitialYPosition();
    vx = planet.getInitialXVelocity();
    vy = planet.getInitialYVelocity();
    r = sqrt(x*x + y*y);
    pi = acos(-1.0);
    FourPi2 = 4.*pi*pi;
    writeTofile(time, x, y, vx/pi, vy/pi);

    while (time < finalTime){
        x += step*vx;
        y += step*vy;
        vx -= step*FourPi2*x/(r*r*r);
        vy -= step*FourPi2*y/(r*r*r);
        r = sqrt(x*x + y*y);
        time += step;
        writeTofile(time, x, y, vx/pi, vy/pi);
    }
    ofile.close();
}

void Solver:: writeTofile(double time_, double x_, double y_, double vx_, double vy_ )
{
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    ofile << setw(15) << setprecision(8) << time_ << ", ";
    ofile << setw(15) << setprecision(8) << x_ << ", ";
    ofile << setw(15) << setprecision(8) << y_ << ", ";
    ofile << setw(15) << setprecision(8) << vx_ << ", ";
    ofile << setw(15) << setprecision(8) << vy_ << endl;
}

