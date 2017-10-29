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

using std::vector;
using namespace std;

class IsingModel
{
private:
    string outfileName;
public:
    IsingModel();
    IsingModel(string fileName_);

    // Function to initialise energy and magnetization
    void initialize(int, double, int **, double&, double&, bool orderingFixed, long& idum);
    // The metropolis algorithm
    void Metropolis(int, long&, int **, double&, double&, double *, double temperature);
    // prints to file the results of the calculations
    void output(int, int, double, double *);
    // inline function for periodic boundary conditions
    inline int periodic(int i, int limit, int add) {
      return (i+limit+add) % (limit);
    }

};

#endif // ISINGMODEL_H
