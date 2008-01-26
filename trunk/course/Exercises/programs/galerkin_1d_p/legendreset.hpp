/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2006-02-20

  Copyright (C) 2007 Université Joseph Fourier Grenoble 1

  This library is free software; you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License as
  published by the Free Software Foundation; either version 2.1 of the
  License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
  02110-1301 USA
*/
/**
   \file legendre.hpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2008-01-24
 */
#ifndef __LegendreSet_H
#define __LegendreSet_H 1

#include <ublas.hpp>
#include <jacobi.hpp>


namespace Life
{
  /**
   * \class LegendreSet
   * \brief LegendreSet polynomial orthonormal basis
   *
   * This class represents the LegendreSet polynomials up to degree \c
   * Degree on a simplex in dimension \c Dim.
   *
   *
   * The legendre polynomials in 1D, the segment \f$[-1;1]\f$ are defined
   * using Jacobi polynomials as follows:
   * \f$ \phi_i(x) = P_i^{0,0}(x) \f$
   *where \f$P_i^{0,0}(x)\f$ is the i-th Jacobi polynomial evaluated at
   * \f$x \in [-1;1]\f$ with weights \f$(0,0)\f$.
   *
   * \ingroup Polynomial
   * @author Christophe Prud'homme
   *
   * @see G.E. Karniadakis and S.J. Sherwin, ''Spectral/hp Element
   * Methods for CFD,'' Oxford University Press, March 1999.
   *
   */
  template<typename T = double>
  class LegendreSet
  {
  public:

    /** @name Typedefs
     */
    //@{

    typedef LegendreSet<T> self_type;
    typedef T value_type;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    LegendreSet( int N, bool normalize = true )
      :
      M_order( N ),
      M_is_normalized( normalize )
    {
    }
    LegendreSet( LegendreSet const & d )
      :
      M_order( d.M_order ),
      M_is_normalized( d.M_is_normalized )
    {
    }

    ~LegendreSet()
    {}

    //@}

    /** @name Operator overloads
     */
    //@{

    self_type const& operator=( self_type const& d )
    {
      if ( this != &d )
	{
	  M_order = d.M_order;
	  M_is_normalized = d.M_is_normalized;
	}
      return *this;
    }

    matrix_type operator()( node_type const& pt ) const
    {
      points_type pts( pt.size(), 1 );
      ublas::column( pts, 0 ) = pt;
      return evaluate( pts );
    }

    matrix_type operator()( points_type const& pts ) const
    {
      return this->evaluate( pts );
    }

    //@}

    /** @name Accessors
     */
    //@{

    /**
     * Number of polynomials in set
     */
    int size() const { return M_order; }

    /**
     * \return the maximum degree of the LegendreSet polynomial to be
     * constructed
     */
    int degree() const { return M_order; }

    /**
     * \return self as a basis
     */
    self_type const& basis() const { return *this; }

    /**
     * \return true if the LegendreSet polynomials are normalized, false
     * otherwise
     */
    bool isNormalized() const { return M_is_normalized; }

    /**
     * \return the \c familyName()
     */
    std::string familyName() const { return "legendre"; }

    //@}

    /** @name  Mutators
     */
    //@{


    //@}

    /** @name  Methods
     */
    //@{



    /**
     * Evaluation at a set of points of the expansion basis in 2D on
     * the triangle
     */
    matrix_type
    evaluate( points_type const& __pts ) const
    {
      matrix_type m ( JacobiBatchEvaluation<value_type>( M_order, 0.0, 0.0, ublas::row(__pts, 0) ) );
      if ( M_is_normalized )
	{
	  for ( int i = 0;i < m.size1(); ++i )
	    ublas::row( m, i ) *= normalization( i );
	}
      return m;
    }

    /**
     * derivation at a set of points of the expansion basis in 2D on
     * the triangle
     */
    template<typename AE>
    vector_matrix_type
    derivate( ublas::matrix_expression<AE> const& __pts ) const
    {
      vector_matrix_type D( 1 );
      D[0].resize( M_order+1, __pts().size2() );
      D[0] = JacobiBatchDerivation<value_type>( M_order, 0.0, 0.0, ublas::row(__pts(),0) );
      if ( M_is_normalized )
	for ( int i = 0; i < M_order+1; ++i )
	  ublas::row( D[0], i ) *= normalization( i );
      return D;
    }

    /**
     * compute the zeros of the Legendre polynomial of degree \p
     * degree
     * \return a matrix wit hthe coordinates of the points
     */
    matrix_type zeros( int degree ) const
    {
      matrix_type z( 1, degree );
      typedef T value_type;
      ublas::vector<T> xr( degree );
      ublas::vector<T> wr( degree );

      // get weights and nodes for LegendreSet polynomials
      gaussjacobi<T, ublas::vector<T> >( degree, wr, xr );

      ublas::row( z, 0 ) = xr;
      return z;
    }
    //@}

  private:
  private:

    value_type normalization( int i ) const
    {
      return (M_is_normalized?std::sqrt( value_type( i ) + 0.5 ) : value_type(1));
    }



  private:
    int M_order;
    bool M_is_normalized;
  };


}
#endif /* __LegendreSet_H */

