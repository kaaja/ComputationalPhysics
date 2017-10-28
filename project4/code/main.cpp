#include <iostream>
#include <string>
#include "isingmodel.h"
//#include "lib.h"

using namespace std;

void read_input(int&, int&, double&, double&, double&);

int main(int argc, char* argv[])
{

  string outfileName;
  long idum;
  int **spin_matrix, n_spins, mcs;
  double w[17], average[5], initial_temp, final_temp, E, M, temp_step;

  // Read in output file, abort if there are too few command-line arguments
  if( argc <= 1 ){
    cout << "Bad Usage: " << argv[0] <<
      " read also output file on same line" << endl;
    exit(1);
  }
  else{
    outfileName=argv[1];
  }
  IsingModel project4b(outfileName);
  //ofile.open(outfilename);
  //    Read in initial values such as size of lattice, temp and cycles
  read_input(n_spins, mcs, initial_temp, final_temp, temp_step);
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
    project4b.initialize(n_spins, temperature, spin_matrix, E, M);
    // start Monte Carlo computation
    for (int cycles = 1; cycles <= mcs; cycles++){
      project4b.Metropolis(n_spins, idum, spin_matrix, E, M, w);
      // update expectation values
      average[0] += E;    average[1] += E*E;
      average[2] += M;    average[3] += M*M; average[4] += fabs(M);
    }
    // print results
    project4b.output(n_spins, mcs, temperature, average);
  }
  free_matrix((void **) spin_matrix); // free memory
  return 0;
}

// read in input data
void read_input(int& n_spins, int& mcs, double& initial_temp,
        double& final_temp, double& temp_step)
{
  cout << "Number of Monte Carlo trials =";
  cin >> mcs;
  cout << "Lattice size or number of spins (x and y equal) =";
  cin >> n_spins;
  cout << "Initial temperature with dimension energy=";
  cin >> initial_temp;
  cout << "Final temperature with dimension energy=";
  cin >> final_temp;
  cout << "Temperature step with dimension energy=";
  cin >> temp_step;
} // end of function read_input
