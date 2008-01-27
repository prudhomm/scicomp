/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2008-01-27

  Copyright (C) 2008 Université Joseph Fourier (Grneoble 1)

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
   \file functionalset.hpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2008-01-27
 */
#ifndef __FunctionalSet_H
#define __FunctionalSet_H 1

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <ublas.hpp>

#include <functional.hpp>

namespace Life
{
/**
 * \class FunctionalSet
 * \brief Set of functionals
 *
 * \ingroup Polynomial
 * @author Christophe Prud'homme
 * @see
 */
template<typename Space>
class FunctionalSet
{
public:


    /** @name Typedefs
     */
    //@{

    typedef Space space_type;
    typedef typename space_type::value_type value_type;


    typedef FunctionalSet<Space> functionalset_type;
    typedef functionalset_type self_type;
    typedef Functional<Space> functional_type;



    typedef std::vector<functional_type> fset_type;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    FunctionalSet()
        :
        _M_space(),
        _M_fset(),
        _M_mat()
    {}

    FunctionalSet( space_type const& s )
        :
        _M_space( s ),
        _M_fset(),
        _M_mat()
    {
    }
    FunctionalSet( space_type const& s, std::vector<functional_type> const& fset )
        :
        _M_space( s ),
        _M_fset( fset ),
        _M_mat( fset.size(), fset[0].coeff().size2() )
    {
        //std::cout << "FunctionalSet: " << fset[0].coeff() <<  "\n";
        this->setFunctionalSet( fset );
    }
    FunctionalSet( FunctionalSet const & fset )
        :
        _M_space( fset._M_space ),
        _M_fset( fset._M_fset ),
        _M_mat( fset._M_mat )
    {}

    ~FunctionalSet()
    {}

    //@}

    /** @name Operator overloads
     */
    //@{

    self_type& operator=( self_type const& fset )
    {
        if ( this != fset )
        {
            _M_space = fset._M_space;
            _M_fset = fset._M_fset;
            _M_mat = fset._M_mat;
        }
        return *this;
    }

    /**
     * \return the i-th functional
     */
    functional_type const& operator()( int i ) const { return _M_fset[i]; }

    /**
     * \return the value of the functional set applied to a polynomial
     */
    matrix_type operator()( space_type const& p ) const
    {
      return ublas::prod( _M_mat, ublas::trans( p.coeff() ) );
    }
    //@}

    /** @name Accessors
     */
    //@{

    /**
     * \return the function space from which the functionals take
     * their values
     */
    space_type const& functionSpace() const { return _M_space; }

    /**
     * This works only if the basis is orthonormal
     *
     * \return the representation of the functional set using basis
     * of the function space.
     */
    matrix_type const& rep() const { return _M_mat; }


    //@}

    /** @name  Mutators
     */
    //@{

    /**
     * set the function space
     */
    void setFunctionSpace( space_type const& __space ) { _M_space = __space; }

    /**
     * set the Functional set
     */
    void setFunctionalSet( std::vector<functional_type> const& fset )
    {
        _M_fset = fset;

	// update matrix associated with functionals applied to the
	// basis of the function space
	_M_mat = ublas::zero_matrix<value_type>( fset.size(), fset[0].coeff().size2() );
	for ( int i = 0; i < fset.size(); ++i )
	  {
	    ublas::row( _M_mat, i ) = ublas::row( fset[i].coeff(), 0 );
	  }
    }

    //@}

    /** @name  Methods
     */
    //@{


    //@}



protected:

private:
    space_type _M_space;
    fset_type _M_fset;
    matrix_type _M_mat;
};
} // Life
#endif /* __FunctionalSet_H */
