## Simulation of the 1D and 2D diffusion equation.

This repository contains c++ files that do the calculations, and a python file for visualizing the outputs from the c++ files.

## Running
The file runProject5.py can be used to run different cases.  

To get an overwiev of how to run runProject4.py, type the following comands in the terminal:  
> python runProject5.py -h  
  
The above command gives you the necesary commandline input to make the program run.  
  
The possible run alternatives are  
    python runproject4.py 5c  
    python runproject4.py 5d  
    python runproject4.py 5f  
    python runproject4.py 5g  
  
where the different cases do the following:  
  - 5c runs 1D simulation for a coarse mesh and a fine mesh, plots and compares results for Forward Euler, Backward Euler, Crank-Nicolson and analytical solution.  
  - 5d runs 1D simulation and plots L2 norms for three different schemes with coarse and fine mesh.  
  - 5f runs 2D simulation for an explicit and an implicit scheme, plots surfaces for both schemes and the analytical solution.  
  - 5g runs 2D simulaton and compares number of iterations as mesh is refined for the implicit Jacobi scheme and the implicit Gauss-Seidel scheme  

If one wants to run the c++ files directly, not visualizing and only get data, type the following in the terminal:  
  
> c++  -fopenmp -o project5.x main.cpp solver.cpp lib.cpp twodimensionaldiffusionsolver.cpp -larmadillo -llapack -lblas  
  
> ./project5.x outfileName dt dx theta endTime scenario numberOfThreads  
  
Of the above command line inputs, the last two need extra mentioning.  
  
"scenario":  
Type "1D" for simulating the one-dimensional case, and 2D for the two-dimensional case.  
If one types "2DJacobiIterations", we can get the number of iterations for both Jacobi scheme and Gauss-Seidel scheme.  
  
"numberOfThreads":  
Sets number of threads one wants to use while running the program, e.g. 2,4,8.  
