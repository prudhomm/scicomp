#include <fstream>
#include <iostream>
#include <iomanip>

#include <polynomialset.hpp>

int main( int argc, char** argv )
{
 int degree = 2;
  if ( argc == 2 )
    degree = std::atoi( argv[1] );

  Life::OrthonormalPolynomialSet<double> opset( degree );
  matrix_type pts( linspace( -1, 1, 30 ) );
  std::cout << "points : " << pts << "\n";
  matrix_type v = opset.evaluate( pts );
  std::cout << "v : " << v << "\n";

  std::ofstream ofs( "legendre.dat" );
  ofs.precision( 10 );
  ofs.setf( std::ios::scientific );
  for( int i = 0; i < v.size2(); ++i )
    {
      ofs << pts(0,i) << " ";
      for( int j = 0; j < v.size1(); ++j )
	ofs << v(j, i )<< " ";
      ofs << "\n";

    }

}
