# Simulation of planatary orbits

This repository contains c++ files that do the calculations, and a python file for visualizing the outputs from
the c++ files.

## Running
The file runProject3.py can be used to run different cases.  

To get an overwiev of how to run runProject3.py, type the following comands in the terminal:
> python runProject3.py -h

The above command gives you the necesary commandline input and an explanation, to make the program run.  
You will then get a list of argument inputs to run different cases: 
    --twoBodyVelocityVerlet  
    --twoBodyFEuler  
    --threeBodyFixedSun  
    --threeBodyFixedSunJupiter10  
    --threeBodyFixedSunJupiter1000  
    --threeBodyMovingSun  
    --threeBodyMovingSunJupiter10  
    --threeBodyMovingSunJupiter1000  
    --solarSystemFixedSun  
    --solarSystemMovingSun  
    --solarSystemMovingSunInnerPlanets  
    --terminalVelocity  
    --terminalVelocityAlternativeForce  
    --perihelion  

An example would be:  
> python runProject3.py --solarSystemMovingSun

If one wants to run the c++ files directly, and not visualizing, type the following in the terminal:  

> c++  -O3 -o project3.x main.cpp planet.cpp solver.cpp planetgeneralrelativityforce.cpp planetalternativeforce.cpp system.cpp

> ./project3.x outfileName, finalTime, N, solverType, initialVy, beta, scenario

Here you can choose "solverType" to be:  
  VelocityVerlet  
  ForwardEuler  
"scenario" can be either one of the command line inputs to runProject3.py  

An example would be:  
> ./project3.x outfileName 10 100000 VelocityVerlet 6.28 3.0 solarSystemMovingSun
