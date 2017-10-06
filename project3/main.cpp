#include "planet.h"
#include "solver.h"

using namespace std;

int main()
{
    int N = 1000;
    double finalTime = 1.0;
    double pi = acos(-1.0);
    Planet earth(1., 1., 0., 0. ,2*pi);
    Solver solution(N, finalTime, "onePlanetSystem.csv");
    solution.addPlanet(earth);
    //cout << earth.getInitialXPosition() << forwardEuler.
    solution.forwardEuler();
    return 0;

}
