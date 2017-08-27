/* Program for 2nd derivative of atan(sqrt(2)).
 Heavily based on MHJ Ch 3.1. */

using namespace std;
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>

void initialize(double *, double *, int *);
void second_derivative(double *, double *, double, double, int);
void output(double *, double *, double, int, char *);

ofstream ofile;

int main(int argc, char *argv[]){
  // Declaration variables
  char *outfilename;
  int number_of_steps;
  double x, initial_step;
  double *computed_derivative, *h_step;
  
  if( argc <= 1){
    cout << "Give a name to the output file after the program name" << endl;
    exit(1);
  }
  else{
    outfilename=argv[1];
  }

  // Read variables from command line
  initialize (&x, &initial_step, &number_of_steps);

  // Pointer/arrays
  h_step = new double[number_of_steps];
  computed_derivative = new double[number_of_steps];

  second_derivative(computed_derivative, h_step, x, initial_step, number_of_steps);

  output(h_step, computed_derivative, x, number_of_steps, outfilename);

  delete [] h_step;
  delete [] computed_derivative;

  return 0;
}

void initialize( double *x, double *initial_step, int *number_of_steps )
{
  printf("Read in from screen:x,  initial step,  and number of steps\n");
  scanf("%lf %lf %d",x,initial_step, number_of_steps);
  return;
}

void second_derivative(double *computed_derivative, double *h_step, double x, double initial_step, int number_of_steps ){
  double h = initial_step;
  for (int i = 0; i < number_of_steps; i++){
    h_step[i] = h;
    computed_derivative[i] = (atan(x+h) - atan(x))/h;
    h = h*0.5;
  }
  return;
}
void output(double *computed_derivative, double *h_step, double x, int number_of_steps, char *outfile_name){
  ofile.open(outfile_name);
  ofile << "log h,  log rel error" << endl;
  ofile << setiosflags(ios::showpoint | ios::uppercase);
  for( int i = 0; i < number_of_steps; i++){
    ofile << setw(15) << setprecision(8) << log10(h_step[i]) << ", ";
    ofile << setw(15) << setprecision(8) << log10(fabs(computed_derivative[i] - 1./3)/(1./3)) << endl;
  }
  ofile.close();
}


