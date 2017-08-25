// Warm-up from MHJ https://github.com/CompPhysics/ComputationalPhysics/blob/master/doc/Projects/2017/WarmUpExercise/warmup.pdf

using namespace std;

# include <iostream>
# include <cmath>
# include <fstream>
# include <iomanip>
#include <cstdlib>

ofstream ofile;

void second_derivative( int , double ,
			double, double *,
			double *, double *);
//void write_file(double der1, double der2, double der3, double der4);
//double f(double x);

void write_to_file( double *, double *, int, char *);


int main(int argc, char* argv[])
{
  //declarations of variables
  double x = sqrt(atof(argv[1]));
  int number_of_steps = atoi(argv[2]);
  double initial_step = atof(argv[3]);
  char *outfilename;
  outfilename = argv[4];
  double *computed_derivative;
  double *h_step;
  double *relative_error;
  

  computed_derivative = new double[number_of_steps];
  h_step = new double[number_of_steps];
  relative_error = new double[number_of_steps];
  second_derivative(number_of_steps, x,
		     initial_step,  h_step,
		    computed_derivative, relative_error );
  write_to_file(relative_error, h_step, number_of_steps, outfilename);
  // cout << h_step[2] << endl;
  // cout << " comp der " << computed_derivative[2] << endl;  
  delete [] h_step;
  delete [] computed_derivative;
  delete [] relative_error;
  return 0;
}

void second_derivative( int number_of_steps, double x,
			double initial_step, double *h_step,
			double *computed_derivative, double *relative_error )
{
  double h;
  h = initial_step;
  int i;
  for (i = 0; i < number_of_steps; i++){
    computed_derivative[i] = (atan(x + h) - atan(x))/h;
    h = h*0.5;
    h_step[i] = h;
    relative_error[i] = log10(fabs(computed_derivative[i] - 1.0/3.0)/(1.0/3.0));
   }

  //cout << " h_step " << *h_step << endl;
  // cout << " comp der " << computed_derivative[0] << endl;
  // cout <<h<< endl;
  return;
}


void write_to_file( double *relative_error, double *step_size, int number_of_steps, char *outfilename)
{
  ofile.open(outfilename);
  ofile << "step size,relative error"<< endl;
  for (int i = 0; i < number_of_steps; i++)
    {
    ofile << setw(15) << setprecision(8) << step_size[i]  << "," << relative_error[i] <<  endl;
    }
    ofile.close();
}



