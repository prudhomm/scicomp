/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>

  Copyright (C) 2008 Christophe Prud'homme
  Copyright (C) 2008 Université Joseph Fourier (Grenoble 1)

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file polynomialset.hpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2008-01-25
 */
#ifndef __PolynomialSet_H
#define __PolynomialSet_H 1

#include <vector>
#include <boost/foreach.hpp>

#include <legendreset.hpp>
#include <polynomial.hpp>


namespace Life
{

  /**
   * \class PolynomialSet
   *  \brief a Set of polynomials
   *
   * This class represents a set of polynomials \f$ {p_i}_{i=1...N} \f$
   * defined in a certain basis given by the template argument. The
   * coefficients of the polynomials in the basis are represented by a
   * matrix whose line represents the polynomials and columns the basis
   * index \f$ C_{i,j} = \mathcal{R}( p_i)_j \f$ where \f$\mathcal{R}\f$
   * is the mapping between the polynomial and its coefficients.
   * We have that the polynomial set is represented as follows:
   * \f[p_i = \sum_j=1^N \mathcal{R}( p_i )_j \phi_j\f]
   *
   *  \ingroup Polynomial
   *  @author Christophe Prud'homme
   *  @see
   */
  template<typename Poly>
  class PolynomialSet
  {
  public:

    /** @name Constants
     */
    //@{

    //@}

    /** @name Typedefs
     */
    //@{
    typedef PolynomialSet<Poly> self_type;
    typedef typename Poly::value_type value_type;
    typedef typename Poly::basis_type basis_type;
    typedef Polynomial<Poly> polynomial_type;
    //@}

    /** @name Constructors, destructor
     */
    //@{

    PolynomialSet()
      :
      _M_basis(),
      _M_coeff(),
      _M_fname( "pset" )
    {
    }
    PolynomialSet( Poly const& p )
      :
      _M_basis( p.basis() ),
      _M_coeff( p.coeff() ),
      _M_fname( p.familyName() )
    {}
    /**
     */
    PolynomialSet( Poly const& p, matrix_type const& c, bool __as_is = false )
      :
      _M_basis( p.basis() ),
      _M_coeff( p.coeff() ),
      _M_fname( p.familyName() )
    {
      setCoefficient( c, __as_is );
    }

    /**
     */
    PolynomialSet( matrix_type const& c, bool __as_is = false )
      :
      _M_basis(),
      _M_coeff( c ),
      _M_fname( "pset" )
    {
      setCoefficient( c, __as_is );
    }


    PolynomialSet( PolynomialSet const & p )
      :
      _M_basis( p._M_basis ),
      _M_coeff( p._M_coeff ),
      _M_fname( p._M_fname )
    {
    }

    virtual ~PolynomialSet()
    {}

    //@}

    /** @name Operator overloads
     */
    //@{

    self_type& operator=( self_type const& pset )
    {
      if ( this != &pset )
	{
	  _M_basis = pset._M_basis;
	  _M_coeff = pset._M_coeff;
	}
      return *this;
    }

    //@}

    /** @name Accessors
     */
    //@{

    /**
     * \return the geometric dimension
     */
    int dim() const { return 1; }

    /**
     * \return the degree of the polynomials in the set
     */
    int degree() const { return _M_basis.degree(); }

    /**
     * \return coefficient of the polynomials of the polynomial set in
     * the basis associated with Poly
     */
    matrix_type const&  coeff() const { return  _M_coeff; }

    /**
     * \return the basis associated with the polynomial set
     */
    basis_type const& basis() const { return _M_basis; }



    /**
     * \return \c true if the polynomial set is zero, \c false otherwise
     */
    bool isZero() const
    {
      return ublas::norm_frobenius( _M_coeff ) < 1e-10;
    }

    /**
     * the \c familyName() identifies the finite element
     * \return the family name of a finite element
     */
    virtual std::string familyName() const  { return _M_fname; }

    /**
     * the name of a finite element is a string composed by:
     *
     * -# a prefix which identifies the family (eg lagrange)
     * -# the dimension
     * -# the order of the finite element
     *
     * \param sep separator between family name, dimension and order (by default it is '.')
     *
     * \return the name of the finite element
     */
    std::string name( std::string sep = "." ) const
    {
      std::ostringstream os;
      os << this->familyName() << sep << dim() << sep << _M_basis.degree();
      return os.str();
    }

    /**
     * \return the number of degrees of freedom
     */
    int nDof() const { return _M_coeff.size1(); }

    //@}

    /** @name  Mutators
     */
    //@{

    /**
     * set the coefficient of the polynomial set in the basis.
     */
    void setCoefficient( matrix_type const& __c, bool __as_is = false )
    {

	  if ( !__as_is )
	    _M_coeff = ublas::prod( __c, _M_coeff );
	  else
	    _M_coeff = __c;

    }

    //@}

    /** @name  Methods
     */
    //@{

    /**
     * Extract the polynomials whose indices are listed in \p list_p
     *
     * \param list_p list of indices of polynomials to extract
     * \return the polynomial set extracted
     */
    PolynomialSet<Poly> polynomials( std::vector<int> const& list_p  ) const
    {
      int dim_p = this->polynomialDimension();
      int new_dim_p = list_p.size();
      matrix_type coeff( list_p.size(), _M_coeff.size2() );
      int j = 0;
      BOOST_FOREACH( int i, list_p )
	{
	  ublas::project( coeff,
			  ublas::range( j, j+1 ),
			  ublas::range( 0, _M_coeff.size2() ) ) =
	    ublas::project( _M_coeff,
			    ublas::range( i, i+1 ),
			    ublas::range( 0, _M_coeff.size2() ) );
	  ++j;
	}
      return PolynomialSet<Poly>( Poly(), coeff, true );
    }

    /**
     * Extract the polynomials up to dimension \p dim
     *
     * \param dim_p polynomial dimension
     * \return the polynomial set extracted
     */
    PolynomialSet<Poly> polynomialsUpToDimension( int dim_p  ) const
    {
      matrix_type coeff( dim_p, _M_coeff.size2() );
      ublas::project( coeff,
		      ublas::range( 0, dim_p ),
		      ublas::range( 0, _M_coeff.size2() ) ) =
	ublas::project( _M_coeff,
			ublas::range( 0, dim_p ),
			ublas::range( 0, _M_coeff.size2() ) );
      return PolynomialSet<Poly>( Poly(), coeff, true );
    }

    /**
     * Extract the polynomials between dimension \p dim_bot up to dimension \p dim_top
     *
     * \param dim_bot polynomial dimension
     * \param dim_top polynomial dimension
     * \return the polynomial set extracted
     */
    PolynomialSet<Poly>
    polynomialsRange( int dim_bot, int dim_top  ) const
    {
      int dim_p = dim_top-dim_bot;
      matrix_type coeff( dim_p, _M_coeff.size2() );

      ublas::project( coeff,
		      ublas::range( 0, dim_top ),
		      ublas::range( 0, _M_coeff.size2() ) ) =
	ublas::project( _M_coeff,
			ublas::range( dim_bot, dim_top ),
			ublas::range( 0, _M_coeff.size2() ) );
      return PolynomialSet<Poly>( Poly(), coeff, true );
    }



    /**
     * evaluate the i-th polynomial at node __pt
     *
     * \warning this function is not efficient at all, the preferred
     * method is to evaluate at a set of points
     */
    template<typename AE>
    ublas::vector<value_type> evaluate( int i, ublas::vector_expression<AE> const& __pt ) const
    {
      return ublas::row( ublas::prod( _M_coeff, _M_basis( __pt ) ),  i );
    }

    /**
     * evaluate all polynomials at node __pt
     * \return a column matrix
     */
    template<typename AE>
    ublas::matrix<value_type> evaluate( ublas::vector_expression<AE> const& __pt ) const
    {
      return ublas::prod( _M_coeff, _M_basis( __pt ) );
    }


    /**
     * evaluate all polynomials of the set at a set of nodes
     *
     * Constructs \f$A_{i,j} = p_i(x_j) = \sum_{k=1}^N \mathcal{R}(p_i)_k \phi_k(x_j)\f$
     *
     * \arg __pts a column oriented matrix contained the node
     * coordinates (in the columns).
     *
     * \return the matrix \f$A\f$
     */
    template<typename AE>
    matrix_type evaluate( ublas::matrix_expression<AE> const& __pts ) const
    {
      matrix_type m ( _M_basis.evaluate( __pts ) );
      return ublas::prod( _M_coeff, m );
    }


    template<typename AE>
    ublas::vector<matrix_type> derivate( ublas::matrix_expression<AE> const& pts ) const
    {
      ublas::vector<matrix_type> der( _M_basis.derivate( pts ) );
      ublas::vector<matrix_type> res( dim() );
      for ( int i = 0;i < dim(); ++i )
	{
	  res[i].resize( _M_coeff.size1(), pts().size2() );
	  ublas::axpy_prod( _M_coeff, der[i], res[i] );
	}
      return res;
    }

    template<typename AE>
    matrix_type derivate( int i, ublas::matrix_expression<AE> const& pts ) const
    {
      ublas::vector<matrix_type> der( _M_basis.derivate( pts ) );
      matrix_type res( _M_coeff.size1(), pts().size2() );
      ublas::axpy_prod( _M_coeff, der[i], res );
      return res;
    }




    //@}

  private:

    basis_type _M_basis;
    matrix_type _M_coeff;
    std::string _M_fname;
  };


template<typename T>
class OrthonormalPolynomialSet
  :
    public PolynomialSet<LegendreSet<T> >
{
  typedef PolynomialSet<LegendreSet<T> > super;
public:


  typedef OrthonormalPolynomialSet<T> self_type;
  typedef self_type component_basis_type;
  typedef T value_type;
  typedef LegendreSet<T> basis_type;

  OrthonormalPolynomialSet( int N )
    :
    super( basis_type( N ) )

  {
    ublas::matrix<value_type> m( ublas::identity_matrix<value_type>( N+1 ) );
    this->setCoefficient( m, true );
  }
  std::string familyName() const { return "legendre"; }

};

} // Life
#endif /* __PolynomialSet_H */

