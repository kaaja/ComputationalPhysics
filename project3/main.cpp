#include "planet.h"
#include "solver.h"
#include "planetmercury.h"
#include "alternativeplanet.h"

using namespace std;

void makeSolarSystem(string outfileName, Solver &solution);
void initialize ( string& outfileName, double& finalTime, int& N, string& solverType, double &initialVy, double &beta, string& scenario, int argc, char** argv);

int main(int argc, char* argv[])
{
    int N;
    double finalTime, initialVy, beta;
    double initialVyJupiter;
    string outfileName, solverType, scenario;

    initialize( outfileName, finalTime, N, solverType, initialVy, beta, scenario, argc, argv);
    Solver solution(N, finalTime);
    double pi = acos(-1.0);

    if (scenario == "twoBody")
    {
        Planet * sun = new Planet(1., 0., 0., 0. , 0., outfileName, "Sun");
        Planet * earth = new Planet(0.000003, 1., 0., 0. , initialVy, outfileName, "Earth");
        solution.addPlanet(*sun);
        solution.addPlanet(*earth);
    }
    else if (scenario == "perihelion")
    {
        Planet * sun_;
        PlanetMercury * mercury_;
        sun_ = new Planet(1., 0., 0., 0. , 0., outfileName, "Sun");
        mercury_ = new PlanetMercury(1.65E-07, 0.3075, 0., 0. , 12.44, outfileName, "Mercury");
        solution.addPlanet(*sun_);
        solution.addPlanet(*mercury_);
    }
    else if (scenario == "perihelionMovingSun")
    {
        PlanetMercury * sun_;
        PlanetMercury * mercury_;
        sun_ = new PlanetMercury(1., 0., 0., 0. , 0., outfileName, "Sun");
        mercury_ = new PlanetMercury(1.65E-07, 0.3075, 0., 0. , 12.44, outfileName, "Mercury");
        solution.addPlanet(*sun_);
        solution.addPlanet(*mercury_);
        solution.changeToCenterOfMassSystem();
    }
    else if (scenario == "alternativeForce")
    {
        Planet * sun_;
        AlternativePlanet * earth_;
        sun_ = new Planet(1., 0., 0., 0. , 0., outfileName, "Sun");
        earth_ = new AlternativePlanet(0.000003, 1, 0., 0. , initialVy, outfileName, "AlternativeEarth", beta);
        solution.addPlanet(*sun_);
        solution.addPlanet(*earth_);
    }
    else if (scenario == "threeBodies")
    {

        Planet * sun = new Planet(1., 0., 0., 0. , 0., outfileName, "Sun");

        Planet * earth = new Planet(0.000003, 1., 0., 0. , initialVy, outfileName, "Earth");

        double velocityEarth = 29.7559;
        double velocityJupiter = 13.0697;
        double velocityJupiterToEarth = velocityJupiter/velocityEarth;
        initialVyJupiter = 2*pi*velocityJupiterToEarth;

        Planet * jupiter = new Planet(9.5e-4, 5.2, 0., 0., initialVyJupiter, outfileName, "Jupiter"); //9.5e-4

        solution.addPlanet(*sun);
        solution.addPlanet(*earth);
        solution.addPlanet(*jupiter);
    }
    else if (scenario == "threeBodiesJupiterTimes10")
    {
        Planet * sun = new Planet(1., 0., 0., 0. , 0., outfileName, "Sun");
        Planet * earth = new Planet(0.000003, 1., 0., 0. , initialVy, outfileName, "Earth");

        double velocityEarth = 29.7559;
        double velocityJupiter = 13.0697;
        double velocityJupiterToEarth = velocityJupiter/velocityEarth;
        initialVyJupiter = 2*pi*velocityJupiterToEarth;

        Planet * jupiter = new Planet(9.5e-4*10, 5.2, 0., 0., initialVyJupiter, outfileName, "Jupiter");

        solution.addPlanet(*sun);
        solution.addPlanet(*earth);
        solution.addPlanet(*jupiter);

    }
    else if (scenario == "threeBodiesJupiterTimes1000")
    {
        Planet * sun = new Planet(1., 0., 0., 0. , 0., outfileName, "Sun");
        Planet * earth = new Planet(0.000003, 1., 0., 0. , initialVy, outfileName, "Earth");

        double velocityEarth = 29.7559;
        double velocityJupiter = 13.0697;
        double velocityJupiterToEarth = velocityJupiter/velocityEarth;
        initialVyJupiter = 2*pi*velocityJupiterToEarth;

        Planet * jupiter = new Planet(9.5e-4*1000, 5.2, 0., 0., initialVyJupiter, outfileName, "Jupiter");

        solution.addPlanet(*sun);
        solution.addPlanet(*earth);
        solution.addPlanet(*jupiter);

    }
    else if (scenario == "threeBodiesMovingSun")
    {
        //makeSolarSystem(outfileName, solution);
        Planet * sun = new Planet(1., 0., 0., 0. , 0., outfileName, "Sun");

        Planet * earth = new Planet(0.000003, 1., 0., 0. , initialVy, outfileName, "Earth");

        double velocityEarth = 29.7559;
        double velocityJupiter = 13.0697;
        double velocityJupiterToEarth = velocityJupiter/velocityEarth;
        initialVyJupiter = 2*pi*velocityJupiterToEarth;

        Planet * jupiter = new Planet(9.5e-4, 5.2, 0., 0., initialVyJupiter, outfileName, "Jupiter"); //9.5e-4

        solution.addPlanet(*sun);
        solution.addPlanet(*earth);
        solution.addPlanet(*jupiter);

        solution.changeToCenterOfMassSystem();
    }
    else if (scenario == "threeBodiesJupiterMassTimes10MovingSun")
    {
        //makeSolarSystem(outfileName, solution);
        Planet * sun = new Planet(1., 0., 0., 0. , 0., outfileName, "Sun");

        Planet * earth = new Planet(0.000003, 1., 0., 0. , initialVy, outfileName, "Earth");

        double velocityEarth = 29.7559;
        double velocityJupiter = 13.0697;
        double velocityJupiterToEarth = velocityJupiter/velocityEarth;
        initialVyJupiter = 2*pi*velocityJupiterToEarth;

        Planet * jupiter = new Planet(9.5e-3, 5.2, 0., 0., initialVyJupiter, outfileName, "Jupiter"); //9.5e-4

        solution.addPlanet(*sun);
        solution.addPlanet(*earth);
        solution.addPlanet(*jupiter);

        solution.changeToCenterOfMassSystem();
    }

    else if (scenario == "threeBodiesJupiterMassTimes1000MovingSun")
    {
        //makeSolarSystem(outfileName, solution);
        Planet * sun = new Planet(1., 0., 0., 0. , 0., outfileName, "Sun");

        Planet * earth = new Planet(0.000003, 1., 0., 0. , initialVy, outfileName, "Earth");

        double velocityEarth = 29.7559;
        double velocityJupiter = 13.0697;
        double velocityJupiterToEarth = velocityJupiter/velocityEarth;
        initialVyJupiter = 2*pi*velocityJupiterToEarth;

        Planet * jupiter = new Planet(9.5e-1, 5.2, 0., 0., initialVyJupiter, outfileName, "Jupiter"); //9.5e-4

        solution.addPlanet(*sun);
        solution.addPlanet(*earth);
        solution.addPlanet(*jupiter);

        solution.changeToCenterOfMassSystem();
    }
    else if (scenario == "solarSystem")
    {
        makeSolarSystem(outfileName, solution);
    }

    else if (scenario == "solarSystemMovingSun" || "solarSystemMovingSunInnerPlanets")
    {
        makeSolarSystem(outfileName, solution);
        solution.changeToCenterOfMassSystem();
    }

    if (solverType == "ForwardEuler")
        solution.forwardEuler();
    else if (solverType == "VelocityVerlet")
        solution.velocityVerlet();
    //solution.velocityVerlet();
    return 0;
}
void initialize (string& outfileName, double& finalTime, int& N, string& solverType, double &initialVy, double &beta, string& scenario, int argc, char** argv)
{
    if( argc<= 1){
      cout << "Insert: outfile-name, number of simulations, amplification factor, start dimension" << endl;
      exit(1);
    }
    else{
      outfileName=argv[1];
    }
    finalTime = atof(argv[2]);
    N = atoi(argv[3]);
    solverType = argv[4];
    initialVy = atof(argv[5]);
    beta = atof(argv[6]);
    scenario = argv[7];
}

void makeSolarSystem(string outfileName, Solver &solution)
{
    string filename1 = "/home/karl/doc/subj/att/fys4150/project3/dataProject3";
    ifstream ifile(filename1);
    string planet; double mass; double x; double y; double vx; double vy; string planetName2;
    if (ifile.is_open()){
        while ( ifile >> planet >> mass >> x >> y >> vx >> vy >> planetName2)
        {
            Planet * planet_;
            planet_ = new Planet(mass, x, y, vx , vy, outfileName, planetName2);
            solution.addPlanet(*planet_);
        }
        ifile.close();
    }
}
