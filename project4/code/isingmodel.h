#ifndef ISINGMODEL_H
#define ISINGMODEL_H

#include<iostream>
#include<new>
#include<vector>
#include<cmath>
#include <fstream>
#include <iomanip>
#include <string>
#include "lib.h"
#include <cmath>
#include <boost/algorithm/string.hpp>
#include <armadillo>

using std::vector;
using namespace std;
using namespace arma;

class IsingModel
{
private:
    string outfileName;
    int temperatureNumber = 0;
    double w[17];

    void energyDifference(double *w, double temperature); // Array with possible deltaE's

    inline int periodic(int i, int limit, int add) {// Inline function for periodic boundary conditions
      return (i+limit+add) % (limit);
    }

public:
    IsingModel();
    IsingModel(string fileName_);

    // Initialise energy and magnetization
    void initialize(int n_spins, int **spin_matrix,
            double& E, double& M, bool orderingFixed, long& idum);

    // Metropolis algorithm
    void Metropolis(int n_spins, long& idum, int **spin_matrix, double& E, double&M, double temperature, int &acceptedMoves,
                                 int myloop_begin, int myloop_end, double *average, colvec &energyArray);

    // Write to file
    void output(int n_spins, int mcs, double temperature, double *average, int acceptedMoves);

    // Write to file timing results
    void outputMPI(int n_spins, int mcs, int numprocs, double TotalTime );

    // Write to file energy array
    void outputEnergyArray(colvec energyArray, int mcs, double temperature);
};

#endif // ISINGMODEL_H
