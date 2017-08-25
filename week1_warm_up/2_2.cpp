using namespace std;

#include <iostream>



int main( int argc, char* argv[])
{
  double S;
  double S2;
  int N = atoi(argv[1]);
  int i;
  double j = 0;
  // cout << " Choose an N:" << endl;
  // cin >> N;
  for( i=1;i<=N;i++ )
    {
      j = i;
      S = S + 1/j;
    }
  cout << "s = " << S << endl;
  for( int i = N; i>0; i--)
  {
      j = i;
      S2 = S2  +  1/j;
  }
  cout << "S2 = " << S2 << endl;
  
}
    
