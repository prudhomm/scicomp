/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2008-01-27

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
   \file polynomial.hpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2008-01-27
 */
#ifndef __Polynomial_H
#define __Polynomial_H 1

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>


namespace Life
{
  namespace ublas = boost::numeric::ublas;

  template<typename > class PolynomialSet;

  /**
   * \class Polynomial
   * \brief polynomial class
   *
   * The polynomial is expressed in the basis from \p Poly. The
   * coefficients of the polynomial in this basis are represented by a
   * matrix whose lines are the polymomial components coefficients (1 if
   * \code is_scalar == true \endcode, \p nDim if \code is_vectorial ==
   * true\endcode and columns are the basis
   *
   * Evaluating the polynomial at a set of points(or just one point) is
   * then simply a matrix-matrix product.
   *
   * \ingroup Polynomial
   * @author Christophe Prud'homme
   * @see
   */
  template<typename Poly>
  class Polynomial
  {
  public:

    /** @name Constants
     */
    //@{

    //@}


    /** @name Typedefs
     */
    //@{

    typedef Polynomial<Poly> self_type;
    typedef typename Poly::value_type value_type;
    typedef typename Poly::basis_type basis_type;
    typedef matrix_type container_type;


    //@}

    /** @name Constructors, destructor
     */
    //@

    /**
     * default constructor
     */
    Polynomial()
      :
      _M_basis(),
      _M_coeff( _M_basis.coeff() )
    {
    }


    /**
     * constructor giving only the underlying basis
     * \param __poly polynomial whose we take the basis
     */
    Polynomial( Poly const& __poly )
      :
      _M_basis( __poly.basis() ),
      _M_coeff( _M_basis.coeff() )
    {
    }

    /**
     * constructor giving the underlying basis and the coefficient of
     * the polynomial in the basis
     *
     * \param __poly polynomial whose we take the basis
     * \param __coeff coefficients of the polynomial in the basis
     */
    Polynomial( Poly const& __poly, container_type const& __coeff, bool __as_is = false )
      :
      _M_basis( __poly.basis() ),
      _M_coeff( _M_basis.coeff() )
    {
      setCoefficient( __coeff, __as_is );
    }

    /**
     * constructor giving the underlying basis and the coefficient of
     * the polynomial in the basis
     *
     * \param __coeff coefficients of the polynomial in the basis
     */
    Polynomial( container_type const& __coeff, bool __as_is = false )
      :
      _M_basis(),
      _M_coeff( _M_basis.coeff() )
    {
      setCoefficient( __coeff, __as_is );
    }

    /**
     * constructor giving the underlying basis and the coefficient of
     * the polynomial in the basis
     *
     * \param __poly polynomial whose we take the basis
     * \param __coeff coefficients of the polynomial in the basis
     */
    template<class AE>
    Polynomial( Poly const& __poly, ublas::matrix_expression<AE> const& __coeff, bool __as_is = false )
      :
      _M_basis( __poly.basis() ),
      _M_coeff( _M_basis.coeff() )
    {
      setCoefficient( __coeff, __as_is );
    }

    Polynomial( Polynomial const & p )
      :
      _M_basis( p._M_basis ),
      _M_coeff( p._M_coeff )
    {}

    ~Polynomial()
    {}

    //@}

    /** @name Operator overloads
     */
    //@{

    self_type const& operator()( self_type const& __p ) const
    {
      if ( this != &__p )
	{
	  _M_basis = __p._M_basis;
	  _M_coeff = __p._M_coeff;
	}
      return *this;
    }


    /**
     * Evaluate polynomial at point \p __x
     *
     * \param __x the coordinate of the point
     * \return the evaluation of the polynomial at point \p __x
     */
    matrix_type
    operator()( node_type const& __x ) const
    {
      return ublas::prod( _M_coeff, _M_basis( __x ) );
    }

    /**
     * Evaluate polynomial at points \p __pts
     *
     * \param __pts the matrix of coordinates of the points
     * \return the evaluation of the polynomial at point \p __x
     */
    matrix_type operator()( points_type const& __pts ) const
    {
      return ublas::prod( _M_coeff, _M_basis( __pts ) );
    }

    //@}

    /** @name Accessors
     */
    //@{


    /**
     * \return \c true if the polynomial set is zero, \c false otherwise
     */
    bool isZero() const
    {
      return ublas::norm_frobenius( _M_coeff ) < 1e-10;

    }
    /**
     * \return the dof
     */
    matrix_type const& coeff() const { return _M_coeff; }

    /**
     * \return the dof
     */
    matrix_type const& coefficients() const { return _M_coeff; }

    /**
     * \return the basis in which the polynomial is expressed
     */
    basis_type const& basis() const { return _M_basis; }

    //@}

    /** @name  Mutators
     */
    //@{

    /**
     * set the coefficient of the polynomial in the basis.
     */
    void setCoefficient( matrix_type const& __c, bool __as_is = false )
    {
      if ( !__as_is )
	{
	  _M_coeff = ublas::prod( __c, _M_coeff );
	}
      else
	_M_coeff = __c;
    }


    //@}

    /** @name  Methods
     */
    //@{

    /**
     * Evaluate polynomial at point \p __x
     *
     * \param __x the coordinate of the point
     * \return the evaluation of the polynomial at point \p __x
     */
    matrix_type
    evaluate( node_type const& __x ) const
    {
      return ublas::prod( _M_coeff, _M_basis( __x ) );
    }

    /**
     * Evaluate polynomial at points \p __pts
     *
     * \param __pts the matrix of coordinates of the points
     * \return the evaluation of the polynomial at point \p __x
     */
    matrix_type evaluate( points_type const& __pts ) const
    {
      return ublas::prod( _M_coeff, _M_basis( __pts ) );
    }

    template<typename AE>
    matrix_type derivate( int i, ublas::matrix_expression<AE> const& pts ) const
    {
      ublas::vector<matrix_type> der( _M_basis.derivate( pts ) );
      matrix_type res( _M_coeff.size1(), pts().size2() );
      ublas::axpy_prod( _M_coeff, der[i], res );
      return res;
    }

    PolynomialSet<Poly> toSet( bool asis = false ) const
    {
      return PolynomialSet<Poly>( Poly(), _M_coeff, asis );
    }
    //@}


  protected:

  private:

    basis_type _M_basis;
    container_type _M_coeff;
  };

}
#endif /* __Polynomial_H */
