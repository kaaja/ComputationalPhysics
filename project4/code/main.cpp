#include <iostream>
#include <string>
#include "isingmodel.h"
#include "mpi.h"

//#include "lib.h"

using namespace std;

void read_input (string& outfileName, int& n_spins, int& mcs, double& initial_temp,
                 double& final_temp, double& temp_step, bool &orderingFixed, int argc, char** argv);

int main(int argc, char* argv[])
{

  string outfileName;
  long idum;
  int **spin_matrix, n_spins, mcs,  my_rank, numprocs;
  double w[17], average[6], initial_temp, final_temp, E, M, temp_step, total_average[6];
  bool orderingFixed;
  colvec energyArray, total_energyArray;

  //  MPI initializations
    MPI_Init (&argc, &argv);
    MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
    if (my_rank == 0 && argc <= 1) {
      cout << "Bad Usage: " << argv[0] <<
        " read output file" << endl;
      exit(1);
    }
    /*
    if (my_rank == 0 && argc > 1) {
      outfilename=argv[1];
      ofile.open(outfilename);
    }
    */

  // Read in initial values such as size of lattice, temp and cycles
  read_input(outfileName, n_spins, mcs, initial_temp, final_temp, temp_step, orderingFixed, argc, argv);

    int no_intervalls = mcs/numprocs;
    int myloop_begin = my_rank*no_intervalls + 1;
    int myloop_end = (my_rank+1)*no_intervalls;
    if ( (my_rank == numprocs-1) &&( myloop_end < mcs) ) myloop_end = mcs;

    // broadcast to all nodes common variables
    MPI_Bcast (&n_spins, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast (&initial_temp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast (&final_temp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast (&temp_step, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  IsingModel project4b(outfileName);
  spin_matrix = (int**) matrix(n_spins, n_spins, sizeof(int));
  energyArray = zeros<colvec>(mcs);
  total_energyArray = zeros<colvec>(mcs);
  //idum = -1; // random starting point
  idum = -1-my_rank;  // random starting point

  double  TimeStart, TimeEnd, TotalTime;
  TimeStart = MPI_Wtime();

  //    initialise energy and magnetization
  E = M = 0.;
  project4b.initialize(n_spins, spin_matrix, E, M, orderingFixed, idum);

  for ( double temperature = initial_temp; temperature <= final_temp+temp_step; temperature+=temp_step){
    // setup array for possible energy changes
    for( int de =-8; de <= 8; de++) w[de+8] = 0;
    for( int de =-8; de <= 8; de+=4) w[de+8] = exp(-de/temperature);
    // initialise array for expectation values
    for( int i = 0; i < 6; i++) average[i] = 0.;
    for( int i = 0; i < 6; i++) total_average[i] = 0.;

    // start Monte Carlo computation
    int acceptedMoves = 0;
    //for (int cycles = 1; cycles <= mcs; cycles++){
    for (int cycles = myloop_begin; cycles <= myloop_end; cycles++){
      project4b.Metropolis(n_spins, idum, spin_matrix, E, M, w, temperature, acceptedMoves);
      // update expectation values
      average[0] += E;    average[1] += E*E;
      average[2] += M;    average[3] += M*M; average[4] += fabs(M); average[5] += pow(fabs(M),2);
      energyArray(cycles -1) = E; // Remember to uncomment
    }

    // Find total average
    for( int i =0; i < 6; i++){
      MPI_Reduce(&average[i], &total_average[i], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    }

    MPI_Reduce(energyArray.memptr(), total_energyArray.memptr(), mcs, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if ( my_rank == 0) {
      //output(n_spins, mcs, temperature, total_average);
      string outfileNameEnergyArray;
      outfileNameEnergyArray = "Temp" + to_string(temperature);
      boost::erase_all(outfileNameEnergyArray, ".");
      boost::erase_all(outfileNameEnergyArray, "0");
      outfileNameEnergyArray = outfileName + outfileNameEnergyArray + "Mcs" + to_string(mcs) + ".csv";
      total_energyArray.save(outfileNameEnergyArray, csv_ascii); // remember to outcomment
      // print results
      project4b.output(n_spins, mcs, temperature, total_average, acceptedMoves);
    }
  }
  free_matrix((void **) spin_matrix); // free memory

  TimeEnd = MPI_Wtime();
  TotalTime = TimeEnd-TimeStart;
  if ( my_rank == 0) {
      cout << "Time = " <<  TotalTime  << " on number of processors: "  << numprocs  << endl;
      project4b.outputMPI(n_spins, mcs, numprocs, TotalTime);
  }

  // End MPI
  MPI_Finalize ();

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
