# Simulation of phase transitions with the 2D Ising model.

This repository contains c++ files that do the calculations, and a python file for visualizing the outputs from the c++ files.

## Running
The file runProject4.py can be used to run different cases.  

To get an overwiev of how to run runProject4.py, type the following comands in the terminal:  
> python runProject3.py -h

The above command gives you the necesary commandline input to make the program run.  
  
The possible run alternatives are  
    python runproject4.py 4b  
    python runproject4.py 4c  
    python runproject4.py 4d  
    python runproject4.py 4e  

If one wants to run the c++ files directly, and not visualizing, type the following in the terminal:  

> mpic++ -O3 -o project4.x main.cpp lib.cpp isingmodel.cpp -larmadillo -llapack -lblas  
  
> mpirun -n numberOfProcessorsForMPI ./project4.x outfileName numberOfSpins numberOfMCCycles InitialTemperature finalTemperature temperatureStep orderingFixed initializeForEachTemperature printEnergyArray  
 
Of the above command line inputs, the three last need extra mentioning.  
  
"orderingFixed":   
Type "orderingFixed" for keeping all initial spin configuration fixed to give minimum energy.  
Type something else than "orderingFixed" to get a random initial spin configuration.  
  
"initializeForEachTemperature":  
Type "initializeForEachTemperature" if the spin configuration is to be initialized again for every new temperature.  
Type somehting else than "initializeForEachTemperature" to make the intial spin configuration for a new temperature equal the last spin configuration from the previous temperature.  
  
"printEnergyArray":  
Type "printEnergyArray" if you run "4d". This is necessary to get the histograms in 4d.  
If you run something else than "4d", and do not want histogram data, type something else than "printEnergyArray".
