## Main files
Project2 contains three files of importance:  
project2.cpp, which runs an eigenvalue solver of your choosing.  
runProject2.py, which runs project2.cpp, reads output from project2.cpp and makes figures from the project2.cpp output results.  
tests-main.cpp, which runs all the tests written in the file tests-functions for the eigenvalue solvers jacobi.cpp, eigenvalueBisection.cpp and lanczos.cpp

## Running 
To run runProject2.py, run the following comands in the terminal:  
> python runProject1.py -h

The above command gives you the necesary commandline input to make the program run.  
An example would be:  
> python runProject2.py 2 lanczosArmadillo

If one wants to run project2.cpp directly, and not using runProject2.py, type the following in the terminal:  
> c++  -O3 -o project2.x project2.cpp jacobi.cpp eigenvalueBisection.cpp lanczos.cpp -L/usr/local/lib -larmadillo -llapack -lblas  
> ./project2.x outfileName numberOfSimulations amplificationFactor N rhoMax maxIterations tolerance solverType electronType omega convergenceLimit

If one wants to run the tests, type the following in the terminal:  
> c++ -o test.x tests-main.cpp test-functions.cpp jacobi.cpp eigenvalueBisection.cpp lanczos.cpp -L/usr/local/lib -larmadillo -llapack -lblas
> ./test.x
