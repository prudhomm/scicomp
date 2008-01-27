#include <fstream>
#include <iostream>
#include <iomanip>

#include <polynomialset.hpp>
#include <functionalset.hpp>

/* Matrix inversion routine.
   Uses lu_factorize and lu_substitute in uBLAS to invert a matrix */
template<class T>
 bool InvertMatrix (const ublas::matrix<T>& input, ublas::matrix<T>& inverse) {
  using namespace boost::numeric::ublas;
  typedef permutation_matrix<std::size_t> pmatrix;
  // create a working copy of the input
  matrix<T> A(input);
 	// create a permutation matrix for the LU-factorization
  pmatrix pm(A.size1());

  // perform LU-factorization
  int res = lu_factorize(A,pm);
  if( res != 0 ) return false;

  // create identity matrix of "inverse"
  inverse.assign(ublas::identity_matrix<T>(A.size1()));

  // backsubstitute to get the inverse
  lu_substitute(A, pm, inverse);

  return true;
}

int main( int argc, char** argv )
{
 int degree = 2;
  if ( argc == 2 )
    degree = std::atoi( argv[1] );

  // construct the polynomial space spanned by the L2 orthonormal polys
  typedef Life::OrthonormalPolynomialSet<double> space_type;
  space_type opset( degree );

  // construct the points
  matrix_type pts( linspace( -1, 1, degree+1 ) );
  std::cout << "points : " << pts << "\n";

  // construct the set of functional
  Life::PointsEvaluation<space_type> fvec( opset, pts );

  Life::FunctionalSet<space_type> fset( opset, fvec );


  matrix_type lag( fset.rep().size1(), fset.rep().size2()  );
  InvertMatrix( fset.rep(), lag );
  lag = ublas::trans( lag );

  Life::PolynomialSet<space_type> lagrangeset( opset, lag, true );

  matrix_type evalpts( linspace( -1, 1, 30 ) );
  matrix_type v( lagrangeset.evaluate( evalpts ) );

  std::ofstream ofs( "lagrange.dat" );
  ofs.precision( 10 );
  ofs.setf( std::ios::scientific );
  for( int i = 0; i < v.size2(); ++i )
    {
      ofs << evalpts(0,i) << " ";
      for( int j = 0; j < v.size1(); ++j )
	ofs << v(j, i )<< " ";
      ofs << "\n";

    }

}
