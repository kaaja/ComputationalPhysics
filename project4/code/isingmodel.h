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
    void energyDifference(double *w, double temperature);
public:
    IsingModel();
    IsingModel(string fileName_);

    // Function to initialise energy and magnetization
    void initialize(int, int **, double&, double&, bool orderingFixed, long& idum);
    // The metropolis algorithm
    void Metropolis(int, long&, int **, double&, double&, double temperature,
                    int &acceptedMoves, int myloop_begin, int myloop_end, double *average, colvec &energyArray);
    // prints to file the results of the calculations
    void output(int, int, double, double *, int acceptedMoves);
    // inline function for periodic boundary conditions
    inline int periodic(int i, int limit, int add) {
      return (i+limit+add) % (limit);
    }
    void outputMPI(int n_spins, int mcs, int numprocs, double TotalTime );
    void outputEnergyArray(colvec energyArray, int mcs, double temperature);

};

#endif // ISINGMODEL_H
