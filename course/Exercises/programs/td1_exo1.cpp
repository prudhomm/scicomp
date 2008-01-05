#include <iostream>
#include <fstream>
#include <blitz/array.h>

int main()
{
  //N =10 ;
  int N = 4;

  // define a matrix A
  blitz::Array<double,2> A(N, N);
  //A = zeros(10 ,10)
  // set to 0
  A = 0;
  std::cout << "[set to 0] A=" << A << "\n";

  // A = ones(10 ,10)
  // set to 0
  A = 1;
  std::cout << "[set to 1] A=" << A << "\n";

  //A = eye(10 ,10)
  // set to identity
  blitz::firstIndex i;
  blitz::secondIndex j;
  A = !(i-j);
  std::cout << "[set to identity] A=" << A << "\n";

  //A ( N /2 ,:) = -1
  A( N/2, blitz::Range::all() ) = -1;
  std::cout << "[A( N /2 ,:) = -1] A=" << A << "\n";

  //A (1:3 , N /2) = -2
  A (blitz::Range(1,3) , N/2) = -2;
  std::cout << "[A (1:3 , N /2) = -2] A=" << A << "\n";

  //A'
  A.transposeSelf(1,0);
  std::cout << "[A'] A=" << A << "\n";
}
