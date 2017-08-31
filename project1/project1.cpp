/* Project1 fys4150 2017. 
    

*/

using namespace std;
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>

void initialize( int *,int *, double *, double *, double *);
//void output(double *, double *, double, int, char *);
void generate_tridiagonal_matrix(int,double, double, double, double **);
void generate_LU(int, double **, double **, double **);
void LU_tridiagonal_solver(double **, double **, double *, double)

ofstream ofile;

int main(int argc, char *argv[]){
  // Declaration variables
  char *outfilename;
  int number_of_simulations, N;
  double a, b, c;
  double ** computed_tridiagonal_matrix, ** LU_Lower, ** LU_Upper;// *h_step;
  
  if( argc <= 1){
    cout << "Give a name to the output file after the program name" << endl;
    exit(1);
  }
  else{
    outfilename=argv[1];
  }

  // Read variables from command line
  initialize(&number_of_simulations, &N, &a, &b, &c);

  // Pointer/arrays
  // h_step = new double[number_of_simulations];
  computed_tridiagonal_matrix = new double*[N];
  LU_Lower = new double*[N];
  LU_Upper = new double*[N];
  for (int i = 0; i < N; i++ ){
    computed_tridiagonal_matrix[i] = new double[N];
    LU_Lower[i] = new double[N];
    LU_Upper[i] = new double[N];
  }
  generate_tridiagonal_matrix(N, a, b, c, computed_tridiagonal_matrix);
  generate_LU(N, computed_tridiagonal_matrix, LU_Lower, LU_Upper);
  

  cout << "Contents of tridiagonal matrix " << endl;
  for (int  m = 0; m < N; m++ ){
    for (int k = 0; k < N; k++){
      cout << LU_Upper[m][k] << " ";
      }
    cout << endl;
  }

  //output(h_step, computed_derivative, x, number_of_steps, outfilename);

  for (int i = 0; i < N; i++ ){
    delete[] computed_tridiagonal_matrix[i];
    delete[] LU_Lower[i];
    delete[] LU_Upper[i];    
  }
  delete [] computed_tridiagonal_matrix;
  delete [] LU_Lower;
  delete [] LU_Upper;
    
  return 0;
}

void initialize( int *number_of_simulations, int *N, double *a, double *b, double *c )
{
  printf("Read in from screen:number of simulations, matrix dimension, lower diagonal, diagonal, upperdiagonal \n");
  scanf(" %d %d %lf %lf %lf",number_of_simulations, N, a, b, c);
  return;
}

void generate_tridiagonal_matrix(int N,double a, double b, double c, double **computed_tridiagonal_matrix ){
  computed_tridiagonal_matrix[0][1] = c;
  computed_tridiagonal_matrix[0][0] = b;
  computed_tridiagonal_matrix[N-1][N-2] = a;
  computed_tridiagonal_matrix[N-1][N-1] = b;
  
  for (int k = 1; k < N-1; k++){
    computed_tridiagonal_matrix[k][k] = b;
    computed_tridiagonal_matrix[k][k-1] = a;
    computed_tridiagonal_matrix[k][k+1] = c;
    }
  return;
}

void generate_LU(int N, double **computed_tridiagonal_matrix, double **LU_Lower, double **LU_Upper){
  for ( int i = 0; i< N; i++ ){
    LU_Lower[i][i] = 1;
  }
  
  LU_Upper[0][0] = computed_tridiagonal_matrix[0][0];
  //U[row-1][row-1] = A[row-1][row-1];
  for (int i = 0; i < N-1; i++ ){
    if ( i != 0){
      LU_Upper[i][i] = computed_tridiagonal_matrix[i][i] -LU_Lower[i][i-1]*LU_Upper[i-1][i];
    }
    LU_Lower[i+1][i] = computed_tridiagonal_matrix[i+1][i]/LU_Upper[i][i];
    LU_Upper[i][i+1] = computed_tridiagonal_matrix[i][i+1];
    }
  LU_Upper[N-1][N-1] = computed_tridiagonal_matrix[N-1][N-1] - LU_Lower[N-1][N-2]*LU_Upper[N-2][N-1];
}
/*
void output(double *computed_derivative, double *h_step, double x, int number_of_steps, char *outfile_name){
  ofile.open(outfile_name);
  ofile << "log_h,log_rel_error" << endl; // DO NOT USE WHITESPACE BETWEEN VAR-NAMES. Pandas does not handle it.
  ofile << setiosflags(ios::showpoint | ios::uppercase);
  for( int i = 0; i < number_of_steps; i++){
    ofile << setw(15) << setprecision(8) << log10(h_step[i]) << ", ";
    ofile << setw(15) << setprecision(8) << log10(fabs(computed_derivative[i] - 1./3)/(1./3)) << endl;
  }
  ofile.close();
  }*/


