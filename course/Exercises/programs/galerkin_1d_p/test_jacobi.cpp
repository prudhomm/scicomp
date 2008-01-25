#include <cassert>

#include <iostream>
#include <iomanip>

#include <jacobi.hpp>

double f( double const& x )
{
  return x*x*x*x;
}

int main(int argc, char** argv)
{
  int degree = 2;
  if ( argc == 2 )
    degree = std::atoi( argv[1] );
  double res = Life::integrate<double>( degree, f );
  //assert( std::abs( res - 2./5.) < 1e-16 );
  std::cout.precision(16);
  std::cout << "int_{-1,1} x*x = "  << res  << "\n";
  std::cout << "exact result is  " << 2./5.  << "\n";
}
