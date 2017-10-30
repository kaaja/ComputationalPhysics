#include <iostream>
#include <string>
#include "isingmodel.h"
//#include "lib.h"

using namespace std;

void read_input (string& outfileName, int& n_spins, int& mcs, double& initial_temp,
                 double& final_temp, double& temp_step, bool &orderingFixed, int argc, char** argv);

int main(int argc, char* argv[])
{

  string outfileName;
  long idum;
  int **spin_matrix, n_spins, mcs;
  double w[17], average[5], initial_temp, final_temp, E, M, temp_step;
  bool orderingFixed;

  // Read in initial values such as size of lattice, temp and cycles
  read_input(outfileName, n_spins, mcs, initial_temp, final_temp, temp_step, orderingFixed, argc, argv);
  IsingModel project4b(outfileName);
  spin_matrix = (int**) matrix(n_spins, n_spins, sizeof(int));
  idum = -1; // random starting point
  for ( double temperature = initial_temp; temperature <= final_temp; temperature+=temp_step){
    //    initialise energy and magnetization
    E = M = 0.;
    // setup array for possible energy changes
    for( int de =-8; de <= 8; de++) w[de+8] = 0;
    for( int de =-8; de <= 8; de+=4) w[de+8] = exp(-de/temperature);
    // initialise array for expectation values
    for( int i = 0; i < 5; i++) average[i] = 0.;
    project4b.initialize(n_spins, temperature, spin_matrix, E, M, orderingFixed, idum);
    // start Monte Carlo computation
    int acceptedMoves;
    for (int cycles = 1; cycles <= mcs; cycles++){
      acceptedMoves = 0;
      project4b.Metropolis(n_spins, idum, spin_matrix, E, M, w, temperature, acceptedMoves);
      // update expectation values
      average[0] += E;    average[1] += E*E;
      average[2] += M;    average[3] += M*M; average[4] += fabs(M);
    }
    // print results
    project4b.output(n_spins, mcs, temperature, average, acceptedMoves);
  }
  free_matrix((void **) spin_matrix); // free memory
  return 0;
}

// read in input data
void read_input (string& outfileName, int& n_spins, int& mcs, double& initial_temp,
                 double& final_temp, double& temp_step, bool &orderingFixed, int argc, char** argv)
{
    if( argc<= 1){
      cout << "Insert: outfile-name, number of simulations, amplification factor, start dimension" << endl;
      exit(1);
    }
    else{
      outfileName=argv[1];
    }
    n_spins = atoi(argv[2]);
    mcs = atoi(argv[3]);
    initial_temp = atof(argv[4]);
    final_temp = atof(argv[5]);
    temp_step = atof(argv[6]);
    string tempInput = argv[7];
    if (tempInput == "orderingFixed"){
        orderingFixed = true;
    }
    else orderingFixed = false;
}
