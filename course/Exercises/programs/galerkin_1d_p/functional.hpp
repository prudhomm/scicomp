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
   \file functional.hpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2005-10-06
 */
#ifndef __Functional_H
#define __Functional_H 1

#include <boost/operators.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>


namespace Life
{
namespace ublas = boost::numeric::ublas;

/**
 * \class Functional
 * \brief represents a linear functional
 *
 * \f$ \ell: \mathbb{P}_N \mapsto \mathbb{R}, p \rightarrow \ell(p) \f$
 *
 * A functional is defined by a polynomial set and a set of
 * coefficients
 *
 * @author Christophe Prud'homme
 * @see
*/
template<typename Space>
class Functional
    :
        public boost::addable<Functional<Space> >
{
    typedef boost::addable<Functional<Space> > super;
public:


    /** @name Typedefs
     */
    //@{

    typedef typename Space::value_type value_type;
    typedef Functional<Space> self_type;
    typedef Space space_type;
    typedef Space polynomialset_type;

    typedef typename space_type::basis_type basis_type;
    typedef typename space_type::polynomial_type polynomial_type;

    static const int nComponents = 1;

    // representation type for the functionals
    typedef ublas::matrix<value_type> rep_type;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    Functional()
        :
        super(),
        _M_p(),
        _M_coeff()
    {}

    Functional( space_type const& P )
        :
        super(),
        _M_p( P ),
        _M_coeff( _M_p.coeff() )

    {

    }

    //template<class AE>
    Functional( space_type const& P,
                matrix_type const& coeff )
        :
        super(),
        _M_p( P ),
        _M_coeff( coeff )
    {
    }

    Functional( Functional const & __f )
        :
        _M_p( __f._M_p ),
        _M_coeff( __f._M_coeff )
    {
    }

    virtual ~Functional()
    {}

    //@}

    /** @name Operator overloads
     */
    //@{

    self_type& operator=( self_type const& __f )
    {
        if ( this != &__f )
        {
            _M_p = __f._M_p;
            _M_coeff = __f._M_coeff;
        }
        return *this;
    }

    /**
     * add to another functional
     * it generates automatically operator+ thanks to addable
     */
    self_type& operator+=(const self_type& __f )
    {
        _M_coeff += __f.M_coeff;
        return *this;
    }

    /**
     * apply the functional to a polynomial
     *
     *
     * \param p polynomial
     * \return matrix resulting from the application of the functional to the polynomial
     */
    virtual matrix_type operator()( polynomial_type const& p ) const
    {
        return ublas::prod( p.coeff(), ublas::trans(_M_coeff) );
    }

    //@}

    /** @name Accessors
     */
    //@{

    /**
     * \return the dimension of the polynomial space
     */
    int domainDimension() const { return _M_coeff.size2(); }

    /**
     * \return the coefficient of the functional in the basis
     * associated with the polynomial space
     */
    rep_type const& coeff() const { return _M_coeff; }

    //@}

    /** @name  Mutators
     */
    //@{

    void setCoefficient( matrix_type const& __coeff )
    {
        _M_coeff = __coeff;
    }


    //@}

    /** @name  Methods
     */
    //@{


    //@}



protected:

private:
    space_type _M_p;
    rep_type _M_coeff;
};
/**
 * \class PointEvaluation
 * \brief generate the functional associated with a point evaluation
 *
 * Given a function space \f$ E \f$, generate the functional \f$ \ell
 * : E \longrightarrow R \f$ associated with the evaluation
 * of the basis functions of the function space at a point in the
 * geometric space.
 *
 * \author Christophe Prud'homme
 */
template<typename Space>
class PointEvaluation
    :
        public Functional<Space>
{
    typedef Functional<Space> super;
public:

    typedef PointEvaluation<Space> self_type;
    typedef typename super::space_type space_type;
    typedef typename space_type::value_type value_type;

    PointEvaluation()
        :
        super()
    {}
    PointEvaluation( space_type const& b, node_type const& __pt )
        : super( b, ublas::trans( b.basis()( __pt ) ) )
    {
    }
};

/**
 * \class PointsEvaluation
 * \brief generate the functionals associated with  point set
 *
 * Given a function space \f$ E \f$ and a set of points in the
 * geometric space \f$ \{x_i\}_{i=1...N} \f$, generate the set of
 * functionals
 * \f{eqnarray*}
 * \ell_i :& E \longrightarrow R, i = 1...N\\
 *         & f \longrightarrow f(x_i)
 * \f}
 * associated with the evaluation of the basis functions of the
 * function space at the set of points in the geometric space.
 *
 * \author Christophe Prud'homme
 */
template<typename Space>
struct PointsEvaluation
    :
        public std::vector<Functional<Space> >
{
    typedef std::vector<Functional<Space> > super;
public:

    typedef PointsEvaluation<Space> self_type;
    typedef Space space_type;


    PointsEvaluation()
        :
        super()
    {}
    PointsEvaluation( space_type const& b, points_type const& __pts )
        : super()
    {
        for ( int c = 0; c < __pts.size2(); ++c )
        {
            //std::cout << "[PointsEvaluation] eval at point " << ublas::column( __pts, c)  << "\n";
          this->push_back( PointEvaluation<Space>( b, ublas::column( __pts, c) ) );
        }

    }
}; // PointsEvaluation



}
#endif /* __Functional_H */
