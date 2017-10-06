#include "solver.h"

Solver:: Solver() { N = finalTime = 0;}
Solver:: Solver(int N_, double finalTime_)
{
    N = N_, finalTime = finalTime_;
    step = finalTime/double(N);
    time = 0.0;

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

    while (time < finalTime){
        x += step*vx;
        y += step*vy;
        vx -= step*FourPi2*x/(r*r*r);
        vy -= step*FourPi2*y/(r*r*r);
        r = sqrt(x*x + y*y);
        time += step;
        cout << "t " << time << " x " << x << " y " << y << " vx " << vx << " vy " << vy  << endl;
    }
}


