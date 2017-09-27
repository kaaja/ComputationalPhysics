#include "jacobi.h"
#include "eigenvalueBisection.h"
#include "time.h"
#include <fstream>
#include <iomanip>


void initialize(string& outfile_name, int& number_of_simulations,int& amplificationFactor, int& N, double& rhoMax,int& maxIterations, double& tolerance, string& solverType, string& interactionRepulsion, double& omega, double& convergenceLimit, int argc, char** argv );
void output_scalars( double computedError, double h, double timeUsed, int N, int counter, double rhoMax, double omega, colvec eigenValues, string convergenceSuccess);
void createTridiagonalMatrix( mat &A, int N, double rhoMax, double rhoMin, double *h, string& interactionRepulsion, double omega);
void calculateError(colvec eigenValues, double *computedError);


ofstream ofile1; // File for scalars

main(int argc, char* argv[]){
  double tolerance, computedError, h, timeUsed, omega, rhoMax, convergenceLimit, lastEigenvalue;
  int N, amplificationFactor, numberOfSimulations, maxIterations, counter;
  string outfileName, outfileNameScalars, solverType, interactionRepulsion, convergenceSuccess;
  mat A, v, eigenvectorMatrixSorted3;
  colvec eigenValues;

  clock_t start, finish;

  initialize(outfileName, numberOfSimulations, amplificationFactor,N, rhoMax, maxIterations, tolerance, solverType, interactionRepulsion, omega, convergenceLimit, argc, argv );

  outfileNameScalars = (outfileName) + string("_scalars")+string(".csv");
  ofile1.open(outfileNameScalars);
  ofile1 << "rhoMax,omega,h,logH,relError,logRelError,timeUsed,logTimeUsed,N,logN,counter,logCounter,lambda1,lambda2,lambda3,convergenceSuccess" << endl;

  // Solving for different matrix dimensions
  eigenValues  = zeros<colvec>(2);
  lastEigenvalue = 5;
  convergenceSuccess = "False";
  for (int simulationNumber = 0; simulationNumber < numberOfSimulations; simulationNumber++){

      //cout << "Rel eigError: " << abs((eigenValues(0) - lastEigenvalue)/lastEigenvalue) << endl;
      lastEigenvalue = eigenValues(0);
      createTridiagonalMatrix(A, N, rhoMax, 0, &h, interactionRepulsion, omega);
      v = eye<mat>(N,N);
      eigenValues  = zeros<colvec>(N);

      start = clock();
      if (solverType == "jacobi")
          jacobi(A, eigenValues, tolerance, maxIterations, N, &counter, v);
      else if (solverType == "armadillo")
          eig_sym(eigenValues, A);
      else if (solverType == "bisection")
          eigenValues = eigenvals3(A, N, tolerance, maxIterations, 0);
      else if (solverType == "bisectionRevised")
          eigenValues = eigenvals3(A, N, tolerance, maxIterations, 1);
      else {
          cout << "choose solvertype " << endl;
          break;
      }

      finish = clock();
      timeUsed = (double)((finish - start)/double(CLOCKS_PER_SEC));

      calculateError(eigenValues, &computedError);

      if (abs((eigenValues(0) - lastEigenvalue)/lastEigenvalue) < convergenceLimit){
          convergenceSuccess = "True";
          if(solverType == "Jacobi")
              output_scalars(computedError, h, timeUsed, N, counter, rhoMax, omega, eigenValues, convergenceSuccess);
          else{
              output_scalars(computedError, h, timeUsed, N, 0, rhoMax, omega, eigenValues, convergenceSuccess);
          }
          break;
      }

      else {
          if(solverType == "jacobi")
              output_scalars(computedError, h, timeUsed, N, counter, rhoMax, omega, eigenValues, convergenceSuccess);
          else{
              output_scalars(computedError, h, timeUsed, N, 0, rhoMax, omega, eigenValues, convergenceSuccess);
          }
      }
      if (simulationNumber < numberOfSimulations -1)
          N *= amplificationFactor;
  }


  outfileName = (outfileName) + string("_eigenVector.csv");
  eigenvectorMatrixSorted3 = get_eigenvecs(A, v, eigenValues, N);
  eigenvectorMatrixSorted3.save(outfileName, csv_ascii);

  ofile1.close();

  return 0;
}

void initialize(string& outfile_name, int& number_of_simulations,int& amplificationFactor, int& N, double& rhoMax,int& maxIterations, double& tolerance, string& solverType, string& interactionRepulsion, double& omega, double& convergenceLimit, int argc, char** argv )
{
    if( argc<= 1){
      cout << "Insert: outfile-name, number of simulations, amplification factor, start dimension" << endl;
      exit(1);
    }
    else{
      outfile_name=argv[1];
    }
    number_of_simulations = atoi(argv[2]);
    amplificationFactor = atoi(argv[3]);
    N = atoi(argv[4]);
    rhoMax = atof(argv[5]);
    maxIterations = atoi(argv[6]);
    tolerance = atof(argv[7]);
    solverType = argv[8];
    interactionRepulsion = argv[9];
    omega = atof(argv[10]);
    convergenceLimit = atof(argv[11]);
}

//ofile1 << "logH,logRelError,timeUsed,logTimeUsed,N,logN,counter,logCounter" << endl;
void output_scalars( double computedError, double h, double timeUsed, int N, int counter, double rhoMax, double omega, colvec eigenValues, string convergenceSuccess){
  ofile1 << setiosflags(ios::showpoint | ios::uppercase);
  ofile1 << setw(15) << setprecision(16) << rhoMax << ", ";
  ofile1 << setw(15) << setprecision(16) << omega << ", ";
  ofile1 << setw(15) << setprecision(16) << h << ", ";
  ofile1 << setw(15) << setprecision(16) << log2(h) << ", ";
  ofile1 << setw(15) << setprecision(16) << computedError << ", ";
  ofile1 << setw(15) << setprecision(16) << log2(computedError) << ", ";
  ofile1 << setw(15) << setprecision(16) << timeUsed << ", ";
  ofile1 << setw(15) << setprecision(16) << log2(timeUsed) << ", ";
  ofile1 << setw(15) << setprecision(16) << N << ", ";
  ofile1 << setw(15) << setprecision(16) << log2(N) << ", ";
  ofile1 << setw(15) << setprecision(16) << counter << ", ";
  ofile1 << setw(15) << setprecision(16) << log2(counter) << ", ";
  ofile1 << setw(15) << setprecision(16) << eigenValues(0)<< ", ";
  ofile1 << setw(15) << setprecision(16) << eigenValues(1) << ", ";
  ofile1 << setw(15) << setprecision(16) << eigenValues(2) << ", ";
  ofile1 << convergenceSuccess << endl;
}
//"rhoMax,omega, h,logH,relError,logRelError,timeUsed,logTimeUsed,N,logN,counter,logCounter,lambda1,lambda2,lambda3"

void createTridiagonalMatrix( mat &A, int N, double rhoMax, double rhoMin, double *h, string& interactionRepulsion, double omega){
    double hTemp = (rhoMax - rhoMin)/N;
    //A.zeros(N,N);
    A = zeros<mat>(N,N);
    double offDiagonal = -1.0/(hTemp*hTemp);
    double diagonalFirstTerm = 2.0/(hTemp*hTemp);

    A(0, 1) =  offDiagonal;
    A(N-1, N-2) = offDiagonal;

    if ( interactionRepulsion == "TwoElectronNoCoulomb" )
        A(0, 0) = diagonalFirstTerm + pow(omega, 2)*hTemp*hTemp;
    else if (interactionRepulsion == "TwoElectronCoulomb" )
        A(0, 0) = diagonalFirstTerm + pow(omega, 2)*hTemp*hTemp + 1./hTemp;
    else
        A(0, 0) = diagonalFirstTerm + hTemp*hTemp;

    for (int row = 1; row < N-1; row++){
        if ( interactionRepulsion == "TwoElectronNoCoulomb" )
            A(row, row) = diagonalFirstTerm + pow(omega, 2)*pow((row+1)*hTemp,2);
        else if (interactionRepulsion == "TwoElectronCoulomb" )
            A(row, row) = diagonalFirstTerm + pow(omega, 2)*pow((row+1)*hTemp,2) + 1./((row+1)*hTemp); // Check this! AM avoiding rho=0 on the first step.
        else
            A(row, row) = diagonalFirstTerm + pow((row+1)*hTemp,2);
      A(row, row - 1) = offDiagonal;
      A(row, row + 1) =  offDiagonal;
      }

    if ( interactionRepulsion == "TwoElectronNoCoulomb" )
         A(N-1, N-1) = diagonalFirstTerm + pow(omega, 2)*pow(rhoMax,2);
    else if (interactionRepulsion == "TwoElectronCoulomb" )
        A(N-1, N-1) = diagonalFirstTerm + pow(omega, 2)*pow(rhoMax,2) + 1./rhoMax;
    else
        A(N-1, N-1) = diagonalFirstTerm + pow(rhoMax,2);
    *h = hTemp;
}

void calculateError(colvec eigenValues, double *computedError){
    // Calculates sup-norm for relative error of 3 lowest eigenvalues
    colvec analyticalSolution  = zeros<colvec>(3);
    analyticalSolution(0) = 3.0;
    analyticalSolution(1) = 7.0;
    analyticalSolution(2) = 11.0;
    *computedError = fabs((analyticalSolution(0) - eigenValues(0))/analyticalSolution(0));;
    double temp;
    for (int i = 1; i < 3; i++) {
        temp = fabs((analyticalSolution(i) - eigenValues(i))/analyticalSolution(i));
        if (temp > *computedError)
            *computedError = temp;
    }
}
