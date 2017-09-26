#include "jacobi.h"

void output_scalars( double computedError, double h, double timeUsed, int N, int counter, double rhoMax, double omega, colvec eigenValues, string convergenceSuccess);
void output_vectors( double *, int, int, string);

ofstream ofile1; // File for scalars
ofstream ofile2; // File for vectors

main(int argc, char* argv[]){
  double tolerance, computedError, h, timeUsed, omega, rhoMax, convergenceLimit, lastEigenvalue;
  int N, amplificationFactor, numberOfSimulations, maxIterations, counter;
  string outfileName, outfileNameScalars, armadillo, interactionRepulsion, convergenceSuccess;
  mat A, v, eigenvectorMatrixSorted3;
  colvec eigenValues;

  clock_t start, finish;

  initialize(outfileName, numberOfSimulations, amplificationFactor,N, rhoMax, maxIterations, tolerance, armadillo, interactionRepulsion, omega, convergenceLimit, argc, argv );

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
      if (armadillo == "false")
          jacobi(A, eigenValues, tolerance, maxIterations, N, &counter, v);
      else
          eig_sym(eigenValues, A);
      finish = clock();
      timeUsed = (double)((finish - start)/double(CLOCKS_PER_SEC));

      calculateError(eigenValues, &computedError);

      if (abs((eigenValues(0) - lastEigenvalue)/lastEigenvalue) < convergenceLimit){
          convergenceSuccess = "True";
          if(armadillo == "false")
              output_scalars(computedError, h, timeUsed, N, counter, rhoMax, omega, eigenValues, convergenceSuccess);
          else{
              output_scalars(computedError, h, timeUsed, N, 0, rhoMax, omega, eigenValues, convergenceSuccess);
          }
          break;
      }

      else {
          if(armadillo == "false")
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
  ofile2.close();

  return 0;
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

void output_vectors( double *output_array, int simulation_number, int N, string outfile_name){
  string filename = outfile_name + to_string(simulation_number+1)+string(".csv");
  ofile2.open(filename);
  /*
  ofile2 << "simulation_number_" << simulation_number+1 << endl;
  ofile2 << setiosflags(ios::showpoint | ios::uppercase);
  for (int i = 0; i<N; i++){
    ofile2 << setw(15) << setprecision(16) << output_array[i] << endl;;
  }
  ofile2 << setw(15) << setprecision(16) << endl;
  */
  ofile2.close();
}
