// Testing of LU-decomposition of tridiagonal system

using namespace std;

#include <iostream>


int main()
{
  int k, m, row = 10 , col = 10;

  double ** A;
  double ** U;
  double ** L;

  // Generating matrix
  A = new double*[row];
  for (int i = 0; i < row; i++ )
    A[i] = new double[col];
  
  A[0][1] = 1;
  A[0][0] = -2;
  A[row-1][row-2] = 1;
  A[row-1][row-1] = -2;
  
  for (k = 1; k < col-1; k++){
    A[k][k] = -2.;
    A[k][k-1] = 1.;
    A[k][k+1] = 1.;
    }

  cout << "Contents of matrix A " << endl;
  for ( m = 0; m < row; m++ ){
    for (k = 0; k < col; k++){
      cout << A[m][k] << " ";
      }
    cout << endl;
  }

  
  // LU decomposition for tridiagonal 
  L = new double*[row];
  U = new double*[row];
  
  for (int i = 0; i < row; i++ ){
    L[i] = new double[col];
    U[i] = new double[col];
  }

  for ( int i = 0; i< row; i++ ){
    L[i][i] = 1;
  }
  
  U[0][0] = A[0][0];
  //U[row-1][row-1] = A[row-1][row-1];
  for (int i = 0; i < row-1; i++ ){
    if ( i != 0){
      U[i][i] = A[i][i] -L[i][i-1]*U[i-1][i];
    }
    L[i+1][i] = A[i+1][i]/U[i][i];
    U[i][i+1] = A[i][i+1];
    }
  U[row-1][row-1] = A[row-1][row-1] - L[row-1][row-2]*U[row-2][row-1];
  
  cout << "Contents of U  " << endl;
  for ( m = 0; m < row; m++ ){
    for (k = 0; k < col; k++){
      cout << U[m][k] << " ";
      }
    cout << endl;
  }

  cout << "Contents of L " << endl;
  for ( m = 0; m < row; m++ ){
    for (k = 0; k < col; k++){
      cout << L[m][k] << " ";
      }
    cout << endl;
  }

  // Freeing memory
  for (int i = 0; i < row; i++ ){
    delete[] A[i];
    delete[] U[i];
    delete[] L[i];
  }
  delete[] A;
  delete[] U;
  delete[] L;
  return 0;
}
