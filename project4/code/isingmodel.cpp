#include "isingmodel.h"

ofstream ofile;


IsingModel::IsingModel(string fileName_)
{
    outfileName = fileName_;
}
// function to initialise energy, spin matrix and magnetization
void IsingModel:: initialize(int n_spins, double temperature, int **spin_matrix,
        double& E, double& M, bool orderingFixed, long& idum)
{
  // setup spin matrix and intial magnetization
  for(int y =0; y < n_spins; y++) {
    for (int x= 0; x < n_spins; x++){
      if(orderingFixed)
        spin_matrix[y][x] = 1; // spin orientation for the ground state
      else
      {
          int spin = (int) (ran1(&idum));
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

void IsingModel:: Metropolis(int n_spins, long& idum, int **spin_matrix, double& E, double&M, double *w, double temperature)
{
  // loop over all spins
  for(int y =0; y < n_spins; y++) {
    for (int x= 0; x < n_spins; x++){
      int ix = (int) (ran1(&idum)*(double)n_spins);
      int iy = (int) (ran1(&idum)*(double)n_spins);
      int deltaE =  2*spin_matrix[iy][ix]*
    (spin_matrix[iy][periodic(ix,n_spins,-1)]+
     spin_matrix[periodic(iy,n_spins,-1)][ix] +
     spin_matrix[iy][periodic(ix,n_spins,1)] +
     spin_matrix[periodic(iy,n_spins,1)][ix]);
      if (ran1(&idum) <= exp(-deltaE/temperature)){//( ran1(&idum) <= w[deltaE+8] ) {
        spin_matrix[iy][ix] *= -1;  // flip one spin and accept new spin config
        M += (double) 2*spin_matrix[iy][ix];
        E += (double) deltaE;
      }
    }
  }
} // end of Metropolis sampling over spins

void IsingModel:: output(int n_spins, int mcs, double temperature, double *average)
{
  double norm = 1/((double) (mcs));  // divided by total number of cycles
  double Eaverage = average[0]*norm;
  double E2average = average[1]*norm;
  double Maverage = average[2]*norm;
  double M2average = average[3]*norm;
  double Mabsaverage = average[4]*norm;
  // all expectation values are per spin, divide by 1/n_spins/n_spins
  double Evariance = (E2average- Eaverage*Eaverage)/n_spins/n_spins;
  double Mvariance = (M2average - Mabsaverage*Mabsaverage)/n_spins/n_spins;
  double Cv = Evariance/(1.*temperature*temperature);
  double chi = Mvariance/(1.*temperature);
  ofile.open(outfileName, ios::out | ios::app);
  ofile << setiosflags(ios::showpoint | ios::uppercase);
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
} // end output function
