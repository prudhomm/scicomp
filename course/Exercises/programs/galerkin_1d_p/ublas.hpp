/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2008-01-21

  Copyright (C) 2008 Université Joseph Fourier

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
#ifndef __ublas_H
#define __ublas_H 1

#include <cmath>

#include <boost/function.hpp>



#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/operation.hpp>

namespace ublas = boost::numeric::ublas;

typedef ublas::vector<double> node_type;
typedef ublas::vector<double> point_type;
typedef ublas::matrix<double,ublas::row_major> points_type;
typedef ublas::matrix<double,ublas::row_major> matrix_type;
typedef ublas::vector<ublas::matrix<double,ublas::row_major> > vector_matrix_type;

/**
 * Construct a set of \p N equidistributed points over \f$[a,b\f$.
 *
 * \param a inf bound of the interval
 * \param b sup bound of the interval
 * \param interior false if include boundary points, false otherwise
 *
 * \return an equidistributed point set over the interval \f$[a;b]\f$
 */
inline
ublas::matrix<double,ublas::row_major>
linspace( double const& a, double const& b, int N, int interior = 0 )
{
    int Npts = N-2*interior;
    ublas::matrix<double,ublas::row_major> v( 1, Npts );
    double h = (b-a)/double(Npts-1);
    double start = a+double(interior)*h;
    for( int i = 0; i < Npts; ++i )
      v(0, i )=start+double(i)*h;
    return v;
}



#endif /* ublas_H */

