using namespace std;
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>

void initialize( char *, int *, int , char *argv[]);

int main(int argc, char *argv[]){
  char filename;
  int  N;
  initialize(&filename, &N, argc, argv);
}

void initialize( char *filename, int *N, int argc, char *argv[] )
{
    if( argc <= 1){
      cout << "Insert: filename, N" << endl;
      exit(1);
    }
    else{
      filename=argv[1];
    }
    *N = atoi(argv[2]);
  return;
}
