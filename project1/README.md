## Dependencies
All files in the directory, excluding the folder 'texPlotsTableData'.

## Main files
Project1 contains two files of importance:  
main.cpp, which does all numerical computations.  
runProject1.py, which runs main.cpp, reads output from main.cpp and plots figures from the main.cpp output results.

## Running the project
To run runProject1.py, run the following comands in the terminal:  
> python runProject1.py -h  

The above command gives you the necesary commandline input to make the program run.  
  
If one wants to run main.cpp directly, and not using runProject.py, type the following in the terminal:
> g++ -o main.x main.cpp lib.cpp
> ./main.x  'outfileName' 'numberOfSimulations' 'amplificationFactor' 'startMeshDimension' 'lowerDiagonal' 'diagonal' 'upperDiagonal'
