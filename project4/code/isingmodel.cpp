#include "isingmodel.h"


ofstream ofile;
ofstream ofile2;


IsingModel::IsingModel(string fileName_)
{
    outfileName = fileName_;
}
// function to initialise energy, spin matrix and magnetization
void IsingModel:: initialize(int n_spins, int **spin_matrix,
        double& E, double& M, bool orderingFixed, long& idum)
{
  // setup spin matrix and intial magnetization
  for(int y =0; y < n_spins; y++) {
    for (int x= 0; x < n_spins; x++){
      if(orderingFixed)
        spin_matrix[y][x] = 1; // spin orientation for the ground state
      else
      {
          int spin = (int) round(ran1(&idum));
          if (spin == 0) spin = -1;
          spin_matrix[y][x] = spin; // spin orientation for the ground state
      }
      M +=  (double) spin_matrix[y][x];
    }
  }
  // setup initial energy
  for(int y =0; y < n_spins; y++) {
    for (int x= 0; x < n_spins; x++){
      E -=  (double) spin_matrix[y][x]*
    (spin_matrix[periodic(y,n_spins,-1)][x] +
     spin_matrix[y][periodic(x,n_spins,-1)]);
    }
  }
}// end function initialise
void IsingModel:: energyDifference(double *w, double temperature)
{
    //for( int de =-8; de <= 8; de+=4) w(de+8) = exp(-de/Temperature);
    for( int de =-8; de <= 8; de++) w[de+8] = 0;
    for( int de =-8; de <= 8; de+=4) w[de+8] = exp(-de/temperature);
}

void IsingModel:: Metropolis(int n_spins, long& idum, int **spin_matrix, double& E, double&M,
                             double *w, double temperature, int &acceptedMoves,
                             int myloop_begin, int myloop_end, double *average, colvec &energyArray)
{
  //colvec w;
  //w = zeros<colvec>(17);
    //double w[17];
  energyDifference(w, temperature);
  // loop over all spins
  for (int cycles = myloop_begin; cycles <= myloop_end; cycles++){
      int allSpins = pow(n_spins, 2);
      for(int y =0; y < allSpins; y++) {
          int ix = (int) (ran1(&idum)*(double)n_spins);
          int iy = (int) (ran1(&idum)*(double)n_spins);
          int deltaE =  2*spin_matrix[iy][ix]*
                (spin_matrix[iy][periodic(ix,n_spins,-1)]+
                spin_matrix[periodic(iy,n_spins,-1)][ix] +
                spin_matrix[iy][periodic(ix,n_spins,1)] +
                spin_matrix[periodic(iy,n_spins,1)][ix]);
         if ( ran1(&idum) <= w[deltaE+8] ) {
            spin_matrix[iy][ix] *= -1;  // flip one spin and accept new spin config
            M += (double) 2*spin_matrix[iy][ix];
            E += (double) deltaE;
            acceptedMoves += 1;
          }
      }
      // update expectation values
      average[0] += E;    average[1] += E*E;
      average[2] += M;    average[3] += M*M; average[4] += fabs(M); average[5] += pow(fabs(M),2);
      //if (printEnergyArray)
        energyArray(cycles -1) = E; // Remember to uncomment
  }

} // end of Metropolis sampling over spins

void IsingModel:: output(int n_spins, int mcs, double temperature, double *average, int acceptedMoves)
{
  string outfileName2;
  outfileName2 = "TempNumber" + to_string(temperatureNumber);
  //boost::erase_all(outfileName2, ".");
  //boost::erase_all(outfileName2, "0");
  outfileName2 = outfileName + outfileName2 + ".csv";
  double norm = 1/((double) (mcs));  // divided by total number of cycles
  double Eaverage = average[0]*norm;
  double E2average = average[1]*norm;
  double Maverage = average[2]*norm;
  double M2average = average[3]*norm;
  double Mabsaverage = average[4]*norm;
  double Mabsaverage2 = average[5]*norm;
  // all expectation values are per spin, divide by 1/n_spins/n_spins
  double Evariance = (E2average- Eaverage*Eaverage)/n_spins/n_spins;
  double Mvariance = (Mabsaverage2 - Mabsaverage*Mabsaverage)/n_spins/n_spins;
  double Cv = Evariance/(1.*temperature*temperature);
  double chi = Mvariance/(1.*temperature);
  ofile.open(outfileName2, ios::out | ios::app);
  ofile << setiosflags(ios::showpoint | ios::uppercase);
  ofile << setw(15) << setprecision(8) << acceptedMoves;
  ofile << setw(15) << setprecision(8) << mcs;
  ofile << setw(15) << setprecision(8) << temperature;
  ofile << setw(15) << setprecision(8) << Eaverage/n_spins/n_spins;
  ofile << setw(15) << setprecision(8) << Evariance;
  ofile << setw(15) << setprecision(8) << Maverage/n_spins/n_spins;
  ofile << setw(15) << setprecision(8) << Mvariance;///temperature;
  ofile << setw(15) << setprecision(8) << Mabsaverage/n_spins/n_spins;
  ofile << setw(15) << setprecision(8) << Cv;
  ofile << setw(15) << setprecision(8) << chi << endl;
  ofile.close();
  temperatureNumber += 1;
} // end output function

void IsingModel:: outputMPI(int n_spins, int mcs, int numprocs, double TotalTime )
{
  string outfileName2;
  outfileName2 = outfileName + "NumProcs" + to_string(numprocs) + ".csv";
  ofile2.open(outfileName2, ios::out | ios::app);
  ofile2 << "n_spins mcs numprocs totalTime" << endl;
  ofile2 << setiosflags(ios::showpoint | ios::uppercase);
  ofile2 << setw(15) << setprecision(8) << n_spins;
  ofile2 << setw(15) << setprecision(8) << mcs;
  ofile2 << setw(15) << setprecision(8) << numprocs;
  ofile2 << setw(15) << setprecision(8) << TotalTime;
  ofile2.close();
}
