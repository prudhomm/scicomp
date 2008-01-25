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
#ifndef __Legendre_H
#define __Legendre_H 1

#include <boost/numeric/ublas/matrix.hpp>


namespace Life
{
/**
 * \class Legendre
 * \brief Legendre polynomial orthonormal basis
 *
 * This class represents the Legendre polynomials up to degree \c
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
template<uint16_type Dim,
         uint16_type Degree,
         typename NormalizationPolicy = Normalized<true>,
         typename T = double,
         template<class> class StoragePolicy = StorageUBlas>
class Legendre : public Basis<LegendreTag<Dim, Degree>, T>
{
    typedef Basis<LegendreTag<Dim, Degree>, T> super;
public:

    typedef LegendreTraits<Dim, Degree, NormalizationPolicy, T, StoragePolicy> traits_type;
    static const uint16_type nDim = Dim;
    static const uint16_type nOrder = Degree;
    static const uint16_type nConvexOrder = nOrder+2;
    static const bool is_normalized = NormalizationPolicy::is_normalized;


    /** @name Typedefs
     */
    //@{

    typedef Legendre<Dim, Degree, NormalizationPolicy, T, StoragePolicy> self_type;

    /*
     * for now can be used only as a basis but we might be interested
     * to have then expressed in other basis like in the moment basis
     */
    typedef self_type basis_type;

    typedef typename traits_type::value_type value_type;

    /*
     * Geometry where the polynomials are defined and constructed
     */
    typedef typename traits_type::convex_type convex_type;
    typedef typename traits_type::reference_convex_type reference_convex_type;

    typedef typename traits_type::diff_pointset_type diff_pointset_type;

    /*
     * storage policy
     */
    typedef typename traits_type::storage_policy storage_policy;
    typedef typename traits_type::matrix_type matrix_type;
    typedef typename traits_type::vector_matrix_type vector_matrix_type;
    typedef typename traits_type::matrix_node_type matrix_node_type;
    typedef typename traits_type::points_type points_type;
    typedef typename traits_type::node_type node_type;


    //@}

    /** @name Constructors, destructor
     */
    //@{

    Legendre()
        :
        super( *this ),
        _M_refconvex(),
        _M_pts( _M_refconvex.makePoints( nDim, 0 ) )
    {
    }
    Legendre( Legendre const & d )
        :
        super( *this ),
        _M_refconvex(),
        _M_pts( d._M_pts )
    {
    }

    ~Legendre()
    {}

    //@}

    /** @name Operator overloads
     */
    //@{

    self_type const& operator=( self_type const& d )
    {
        if ( this != &d )
            {
                _M_pts = d._M_pts;
            }
        return *this;
    }

    //ublas::matrix_column<matrix_type const> operator()( node_type const& pt ) const
    matrix_type operator()( node_type const& pt ) const
    {
        points_type pts( pt.size(), 1 );
        ublas::column( pts, 0 ) = pt;
        return evaluate( pts );
    }

    matrix_type operator()( points_type const& pts ) const
    {
        return evaluate( pts );
    }

    //@}

    /** @name Accessors
     */
    //@{

    /**
     * Number of polynomials in set
     */
    size_type size() const { return convex_type::polyDims( nOrder ); }

    /**
     * \return the maximum degree of the Legendre polynomial to be
     * constructed
     */
    uint16_type degree() const { return nOrder; }

    /**
     * \return self as a basis
     */
    self_type const& basis() const { return *this; }

    /**
     * \return true if the Legendre polynomials are normalized, false
     * otherwise
     */
    bool isNormalized() const { return is_normalized; }

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
     * Legendre polynomials is an orthonormal basis, the coefficients
     * of the polynomials of the basis are the canonical vectors and
     * represented by the identity matrix (lines are polynomials and
     * columns are the polynomial basis )
     *
     * This function is correct only if we use the Legendre polynomials
     * as a basis
     */
    matrix_type coeff() const
    {
        return ublas::identity_matrix<value_type>( reference_convex_type::polyDims( nOrder ), _M_pts.size2() );
    }


    /**
     * evaluate the Legendre polynomials at a set of points \p __pts
     *
     * \arg __x is a set of points
     */
    static matrix_type evaluate( points_type const& __pts )
    {
        return evaluate( __pts, mpl::int_<nDim>() );
    }

    template<typename AE>
    static vector_matrix_type derivate( ublas::matrix_expression<AE>  const& __pts )
    {
        return derivate( __pts, mpl::int_<nDim>() );
    }


    //@}

private:
private:

    static value_type normalization( int i )
    {
      return (is_normalized?math::sqrt( value_type( i ) + 0.5 ) : value_type(1));
    }

    /**
     * Evaluation at a set of points of the expansion basis in 2D on
     * the triangle
     */
    static matrix_type
    evaluate( points_type const& __pts, mpl::int_<1> )
    {
        matrix_type m ( JacobiBatchEvaluation<nOrder,value_type>( 0.0, 0.0, ublas::row(__pts, 0) ) );
        if ( is_normalized )
            {
                for ( uint16_type i = 0;i < m.size1(); ++i )
                    ublas::row( m, i ) *= normalization( i );
            }
        return m;
    }

    /**
     * derivation at a set of points of the expansion basis in 2D on
     * the triangle
     */
    template<typename AE>
    static vector_matrix_type
    derivate( ublas::matrix_expression<AE> const& __pts, mpl::int_<1> )
    {
        vector_matrix_type D( 1 );
        D[0].resize( nOrder+1, __pts().size2() );
        D[0] = JacobiBatchDerivation<nOrder,value_type>( 0.0, 0.0, ublas::row(__pts(),0) );
        if ( is_normalized )
            for ( uint16_type i = 0; i < nOrder+1; ++i )
                ublas::row( D[0], i ) *= normalization( i );
        return D;
    }
};


}
#endif /* __Legendre_H */

