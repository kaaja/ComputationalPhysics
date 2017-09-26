// 1st try on Jacobi algorithm for finding eigenvalues

#include <iostream>
#include <armadillo>
#include <cmath>
#include <fstream>
#include <iomanip>
#include "time.h"

//#define CATCH_CONFIG_RUNNER
//#include "catch.hpp"

using namespace std;
using namespace arma;

void initialize(string& outfile_name, int& number_of_simulations,int& amplificationFactor, int& N, double& rhoMax,int& maxIterations, double& tolerance, string& armadillo, string& interactionRepulsion, double& omega, double& convergenceLimit, int argc, char** argv );
void jacobi(mat &A, colvec &eigenValues, double tolerance, int maxIterations, int N, int *counter, mat& v);
mat get_eigenvecs(mat a, mat v, colvec eigenValues, int N);
double findMaxNonDiagonalElement(mat &A, int *k, int *l, int N);
void rotate(mat &aMatrix, int k, int l, int N, mat& v);
void createEigenvalueVector( mat A, colvec &eigenValues, int N );
void createTridiagonalMatrix( mat &A, int N, double rhoMax, double rhoMin, double *h, string& interactionRepulsion, double omega);
void output_scalars( double computedError, double h, double timeUsed, int N, int counter, double rhoMax, double omega, colvec eigenValues, string convergenceSuccess);
void output_vectors( double *, int, int, string);
void calculateError(colvec eigenValues, double *computedError);



void initialize(string& outfile_name, int& number_of_simulations,int& amplificationFactor, int& N, double& rhoMax,int& maxIterations, double& tolerance, string& armadillo, string& interactionRepulsion, double& omega, double& convergenceLimit, int argc, char** argv )
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
    armadillo = argv[8];
    interactionRepulsion = argv[9];
    omega = atof(argv[10]);
    convergenceLimit = atof(argv[11]);
}


void jacobi(mat &A, colvec &eigenValues, double tolerance, int maxIterations, int N, int *counter, mat& v){
    int k,l;
    *counter = 0;

    double aMaxNonDiagonal = findMaxNonDiagonalElement(A, &k, &l, N);

    while ( fabs(aMaxNonDiagonal) > tolerance && *counter < maxIterations ){
        aMaxNonDiagonal  = findMaxNonDiagonalElement(A, &k, &l, N);
        rotate(A, k, l, N, v);
        *counter += 1;
    }
    createEigenvalueVector( A, eigenValues,  N);
}

mat get_eigenvecs(mat a, mat v, colvec eigenValues, int N){
    // Based on Morten's example
    mat vecs(N,3);
    for(int i=0;i<3;i++){
        for(int j=0;j<N;j++){
            if(a(j,j)==eigenValues(i)){
                for(int k=0;k<N;k++){
                      vecs(k,i)=v(k,j);
                }
             }
         }
    }
    return vecs;
}

void rotate(mat &aMatrix, int k, int l, int N, mat& v){
    double c, s;
    if ( aMatrix(k,l) != 0.0 ) {
        double tau, t;
        tau = (aMatrix(l,l) - aMatrix(k,k))/(2.*aMatrix(k,l));
        if ( tau >= 0 ) {
          t = 1.0/(tau + sqrt(1.0 + tau*tau));
        }
        else {
          t = -1.0/(-tau +sqrt(1.0 + tau*tau));
        }

        c = 1/sqrt(1+t*t);
        s = c*t;
    }
    else {
        c = 1.0;
        s = 0.0;
    }

    double aIk, aIl, aKk, aLl, vik, vil;
    for (int i = 0; i < N; i++){
        if ( i != k && i != l ) {
                aIk = aMatrix(i, k);
                aIl = aMatrix(i, l);
                aMatrix(i,k) = aIk*c - aIl*s;
                aMatrix(k, i) = aMatrix(i, k);
                aMatrix(i, l) = aIl*c + aIk*s;
                aMatrix(l, i) = aMatrix(i, l);
        }
        vik=v(i,k);
        vil=v(i,l);
        v(i,k)=c*vik-s*vil;
        v(i,l)=c*vil+s*vik;
    }
    aKk = aMatrix(k,k);
    aLl = aMatrix(l, l);
    aMatrix(k, k) = aKk*c*c - 2.*aMatrix(k, l)*c*s + aLl*s*s ;
    aMatrix(l, l) = aLl*c*c + 2.*aMatrix(k, l)*c*s + aKk*s*s;
    aMatrix(k,l) = 0.0;
    aMatrix(l,k) = 0.0;
}

double findMaxNonDiagonalElement(mat &A, int *k, int *l, int N){
    double aMaxTemp = 0.0;
    for (int row = 0; row < N; row++){
         for (int col = row + 1; col < N; col++){
             if( row != col && fabs(A(row, col)) >  aMaxTemp ) {
                 *k = row;
                 *l = col;
                 aMaxTemp = fabs(A(row, col));
             }
         }
    }
    return aMaxTemp;
}

void createEigenvalueVector( mat A, colvec &eigenValues, int N ){
    for (int row = 0; row < N; row++){
        for (int col= 0; col < N; col++){
            if ( col == row )
                eigenValues(col) = A(row, col);
        }
    }
    eigenValues = sort(eigenValues);
}

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


