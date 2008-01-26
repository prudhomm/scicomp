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
#include <legendreset.hpp>


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

    /**
     * \brief extract the i-th component of a vectorial polynomial set
     *
     * \return the i-th component of the polynomial set
     */
    component_type operator[]( uint16_type i ) const
    {
      BOOST_STATIC_ASSERT( is_vectorial );
      const int nrows = _M_coeff.size1()/nComponents;
      const int ncols = _M_coeff.size2();
      return component_type( Poly(), ublas::project( _M_coeff,
						     ublas::slice( nrows*i+i, nComponents, nrows/nComponents ),
						     ublas::slice( 0, 1, ncols ) ), true );
    }

    //@}

    /** @name Accessors
     */
    //@{

    /**
     * \return the degree of the polynomials in the set
     */
    uint16_type degree() const { return nOrder; }

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
     * \return true if the polynomial set is scalar, false otherwise
     */
    static bool isScalar() { return is_scalar; }

    /**
     * \return true if the polynomial set is vectorial, false otherwise
     */
    static bool isVectorial() { return is_vectorial; }

    /**
     *
     */
    static uint16_type numberOfComponents() { return nComponents; }

    /**
     * \return the polynomial dimension
     */
    size_type polynomialDimension() const { return _M_coeff.size1()/nComponents; }

    /**
     * \return the polynomial dimension per component
     */
    size_type polynomialDimensionPerComponent() const { return _M_coeff.size2();}

    /**
     * \return \c true if the polynomial set is zero, \c false otherwise
     */
    bool isZero() const
    {
      return ublas::norm_frobenius( _M_coeff ) < Life::type_traits<value_type>::epsilon();
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
      os << this->familyName() << sep << nDim << sep << nOrder;
      return os.str();
    }

    //@}

    /** @name  Mutators
     */
    //@{

    /**
     * set the coefficient of the polynomial set in the basis.
     */
    void setCoefficient( matrix_type const& __c, bool __as_is = false )
    {
      if ( is_scalar )
	{
	  if ( !__as_is )
	    _M_coeff = ublas::prod( __c, _M_coeff );
	  else
	    _M_coeff = __c;
	}
      else
	{
	  if ( !__as_is )
	    {
	      _M_coeff = ublas::prod( __c, polyset_type::toMatrix( _M_coeff ) );
	      _M_coeff = polyset_type::toType( _M_coeff );
	    }
	  else
	    _M_coeff = __c;
	}
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
    PolynomialSet<Poly, PolySetType> polynomials( std::vector<uint16_type> const& list_p  ) const
    {
      size_type dim_p = this->polynomialDimension();
      size_type new_dim_p = nComponents*list_p.size();
      matrix_type coeff( nComponents*nComponents*list_p.size(), _M_coeff.size2() );
      int j = 0;
      BOOST_FOREACH( uint16_type i, list_p )
	{
	  for ( int c = 0; c < nComponents; ++c )
	    {
	      ublas::project( coeff,
			      ublas::range( c*new_dim_p+nComponents*j, c*dim_p+nComponents*j+nComponents ),
			      ublas::range( 0, _M_coeff.size2() ) ) =
		ublas::project( _M_coeff,
				ublas::range( c*dim_p+nComponents*i, c*dim_p+nComponents*i+nComponents ),
				ublas::range( 0, _M_coeff.size2() ) );
	    }
	  ++j;
	}
      return PolynomialSet<Poly, PolySetType>( Poly(), coeff, true );
    }

    /**
     * Extract the polynomials up to dimension \p dim
     *
     * \param dim_p polynomial dimension
     * \return the polynomial set extracted
     */
    PolynomialSet<Poly, PolySetType> polynomialsUpToDimension( int dim_p  ) const
    {
      matrix_type coeff( nComponents*nComponents*dim_p, _M_coeff.size2() );
      for ( int c = 0; c < nComponents; ++c )
	{
	  size_type nc = c*this->polynomialDimension();
	  ublas::project( coeff,
			  ublas::range( c*nComponents*dim_p, (c+1)*nComponents*dim_p ),
			  ublas::range( 0, _M_coeff.size2() ) ) =
	    ublas::project( _M_coeff,
			    ublas::range( nc, nc+nComponents*dim_p ),
			    ublas::range( 0, _M_coeff.size2() ) );
	}
      return PolynomialSet<Poly, PolySetType>( Poly(), coeff, true );
    }

    /**
     * Extract the polynomials between dimension \p dim_bot up to dimension \p dim_top
     *
     * \param dim_bot polynomial dimension
     * \param dim_top polynomial dimension
     * \return the polynomial set extracted
     */
    PolynomialSet<Poly, PolySetType>
    polynomialsRange( uint16_type dim_bot, uint16_type dim_top  ) const
    {
      uint16_type dim_p = dim_top-dim_bot;
      matrix_type coeff( nComponents*nComponents*dim_p,
			 _M_coeff.size2() );

      for ( int c = 0; c < nComponents; ++c )
	{
	  size_type nc = c*this->polynomialDimension();
	  ublas::project( coeff,
			  ublas::range( c*nComponents*dim_bot, (c+1)*nComponents*dim_top ),
			  ublas::range( 0, _M_coeff.size2() ) ) =
	    ublas::project( _M_coeff,
			    ublas::range( nc+dim_bot, nc+nComponents*dim_top ),
			    ublas::range( 0, _M_coeff.size2() ) );
	}
      return PolynomialSet<Poly, PolySetType>( Poly(), coeff, true );
    }

    /**
     * Extract the \p i -th polynomial
     *
     * \param i index of the polynomial to extract
     * \return the polynomial extracted
     */
    Polynomial<Poly, PolySetType> polynomial( uint16_type i  ) const
    {
      size_type dim_p = this->polynomialDimension();
      matrix_type coeff( nComponents*nComponents, _M_coeff.size2() );
      for ( int c = 0; c < nComponents; ++c )
	{
	  ublas::project( coeff,
			  ublas::range( c*nComponents, c*nComponents+nComponents ),
			  ublas::range( 0, _M_coeff.size2() ) ) =
	    ublas::project( _M_coeff,
			    ublas::range( c*dim_p+nComponents*i, c*dim_p+nComponents*i+nComponents ),
			    ublas::range( 0, _M_coeff.size2() ) );
	}
      return Polynomial<Poly, PolySetType> ( Poly(), coeff, true );
    }


    /**
     * evaluate the i-th polynomial at node __pt
     *
     * \warning this function is not efficient at all, the preferred
     * method is to evaluate at a set of points
     */
    template<typename AE>
    ublas::vector<value_type> evaluate( uint16_type i, ublas::vector_expression<AE> const& __pt ) const
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
      LIFE_ASSERT( _M_coeff.size2() == m.size1() )(_M_coeff.size2())(m.size1() ).error("invalid size");
      return ublas::prod( _M_coeff, m );
    }

    /**
     * Derivate with respect to the \f$\ell\f$ coordinates at the
     * nodes where the polynomials basis have been constructed.
     *
     * We construct the matrix \f{eqnarray*} A_{i,j} &= \frac{\partial
     * p_i(x_j)}{\partial x_\ell}\\ &= \sum_{k=1}^N \mathcal{R}(p_i)_k
     * \frac{\partial \phi_k(x_j)}{\partial x_\ell} \f}
     *
     * \arg l the derivation index \f$\ell\f$
     * \arg __pts a column oriented matrix contained the node
     * coordinates (in the columns).
     *
     * \return the matrix \f$A\f$
     */

    /**
     * \brief derivatives of Dubiner polynomials
     * the derivatives are computed at the nodes of the lattice
     *
     * \arg i index of the derivative (0 : x, 1 : y, 2 : z )
     */
    matrix_type const& d( uint16_type i ) const
    {
      return _M_basis.d(i);
    }

    matrix_type d( uint16_type i, uint16_type j ) const
    {
      return ublas::prod( _M_basis.d(i), _M_basis.d(j) );
    }

    /**
     * \brief Derivate with respect to the l-th direction.
     *
     * \return the polynomial set associated with the derivation in the l-direction
     */
    self_type derivate( uint16_type l ) const
    {
      return self_type( Poly(), ublas::prod(  _M_coeff, _M_basis.d( l ) ), true );
    }

    template<typename AE>
    ublas::vector<matrix_type> derivate( ublas::matrix_expression<AE> const& pts ) const
    {
      ublas::vector<matrix_type> der( _M_basis.derivate( pts ) );
      ublas::vector<matrix_type> res( nDim );
      for ( uint16_type i = 0;i < nDim; ++i )
	{
	  res[i].resize( _M_coeff.size1(), pts().size2() );
	  ublas::axpy_prod( _M_coeff, der[i], res[i] );
	}
      return res;
    }

    template<typename AE>
    matrix_type derivate( uint16_type i, ublas::matrix_expression<AE> const& pts ) const
    {
      ublas::vector<matrix_type> der( _M_basis.derivate( pts ) );
      matrix_type res( _M_coeff.size1(), pts().size2() );
      ublas::axpy_prod( _M_coeff, der[i], res );
      return res;
    }

    template<typename AE>
    matrix_type derivate( uint16_type i, uint16_type j, ublas::matrix_expression<AE> const& pts ) const
    {
      //std::cout << "[derivate2] _M_coeff = " << _M_coeff << "\n";
      matrix_type eval( _M_basis.evaluate( pts ) );
      //matrix_type res( _M_coeff.size1(), pts().size2() );
      //ublas::axpy_prod( _M_coeff, der[i], res );
      matrix_type p1 = ublas::prod( _M_coeff, _M_basis.d(i) );
      matrix_type p2 = ublas::prod( p1, _M_basis.d(j) );
      return ublas::prod( p2, eval );
    }
    /**
     * Gradient of the polynomial set
     *
     * Computes the gradient of the polynomial set.
     */
    gradient_polynomialset_type
    gradient() const
    {
      const int n1 = _M_coeff.size1();
      const int n2 = _M_coeff.size2();
      ublas::matrix<value_type> c ( nDim*nDim*n1, n2 );
      c.clear();
      for ( int i = 0; i <nDim; ++i )
	{
	  ublas::project( c,
			  ublas::slice( nDim*n1*i+i, nDim, n1 ),
			  ublas::slice( 0, 1, n2 ) )  = ublas::prod( _M_coeff, _M_basis.d( i ) );
	}
      return gradient_polynomialset_type( c, true );
    }



    /**
     * \return the number of degrees of freedom
     */
    uint16_type nbDof() const { return _M_coeff.size1(); }


    //@}

  private:

    basis_type _M_basis;
    matrix_type _M_coeff;
    std::string _M_fname;
  };

  template<typename Poly,template<uint16_type> class PolySetType> const bool PolynomialSet<Poly,PolySetType>::is_scalar;
  template<typename Poly,template<uint16_type> class PolySetType> const bool PolynomialSet<Poly,PolySetType>::is_vectorial;
  template<typename Poly,template<uint16_type> class PolySetType> const bool PolynomialSet<Poly,PolySetType>::is_tensor2;
  template<typename Poly,template<uint16_type> class PolySetType> const uint16_type PolynomialSet<Poly,PolySetType>::nComponents;
  template<typename Poly,template<uint16_type> class PolySetType> const uint16_type PolynomialSet<Poly,PolySetType>::nComponents1;
  template<typename Poly,template<uint16_type> class PolySetType> const uint16_type PolynomialSet<Poly,PolySetType>::nComponents2;

  /**
   * \class OrthonormalPolynomialSet
   * \brief a set of orthonormal polynomials over a convex
   *
   * On the simplicies we use the Dubiner basis
   *
   */
  template<uint16_type Dim,
	   uint16_type Order,
	   template<uint16_type> class PolySetType = Scalar,
	   typename T = double,
	   template<uint16_type,uint16_type,uint16_type> class Convex = Simplex>
  class OrthonormalPolynomialSet
  {};

template<uint16_type Dim,
         uint16_type Order,
         template<uint16_type> class PolySetType,
         typename T>
class OrthonormalPolynomialSet<Dim, Order, PolySetType, T, Simplex>
  :
    public PolynomialSet<Dubiner<Dim, Order, Normalized<true>, T, StorageUBlas>, PolySetType >
{
  typedef PolynomialSet<Dubiner<Dim, Order, Normalized<true>, T, StorageUBlas>, PolySetType > super;
public:

  static const uint16_type nDim = Dim;
  static const uint16_type nOrder = Order;

  typedef OrthonormalPolynomialSet<Dim, Order, PolySetType, T, Simplex> self_type;
  typedef self_type component_basis_type;

  typedef typename super::polyset_type polyset_type;
  static const bool is_tensor2 = polyset_type::is_tensor2;
  static const bool is_vectorial = polyset_type::is_vectorial;
  static const bool is_scalar = polyset_type::is_scalar;
  static const bool is_continuous = false;
  static const bool is_modal = true;
  static const uint16_type nComponents = polyset_type::nComponents;
  typedef typename super::component_type component_type;

  typedef T value_type;
  typedef Dubiner<Dim, Order, Normalized<true>, T, StorageUBlas> basis_type;
  typedef Simplex<Dim, Order, Dim> convex_type;
  template<int O>
  struct convex
  {
    typedef Simplex<Dim, O, Dim> type;
  };
  typedef Reference<convex_type, nDim, nOrder, nDim, value_type> reference_convex_type;

  typedef typename super::polynomial_type polynomial_type;

  //!< Number of degrees of freedom per vertex
  static const uint16_type nDofPerVertex = convex_type::nbPtsPerVertex;
  //!< Number of degrees  of freedom per edge
  static const uint16_type nDofPerEdge = convex_type::nbPtsPerEdge;
  //!< Number of degrees  of freedom per face
  static const uint16_type nDofPerFace = convex_type::nbPtsPerFace;

  //!< Number of degrees  of freedom per volume
  static const uint16_type nDofPerVolume = convex_type::nbPtsPerVolume;

  static const uint16_type nLocalDof = convex_type::numPoints;

  static const uint16_type nDof = nLocalDof;
  static const uint16_type nNodes = nDof;
  static const uint16_type nDofGrad = super::nDim*nDof;
  static const uint16_type nDofHess = super::nDim*super::nDim*nDof;
  typedef typename matrix_node<value_type>::type points_type;

  OrthonormalPolynomialSet()
    :
    super( basis_type() )

  {
    ublas::matrix<value_type> m( ublas::identity_matrix<value_type>( nComponents*convex_type::polyDims( nOrder ) ) );
    if ( is_tensor2 )
      std::cout << "[orthonormalpolynomialset] m = " << m << "\n";
    if ( !(ublas::norm_frobenius( polyset_type::toMatrix( polyset_type::toType( m ) ) -
				  m ) < 1e-10 ) )
      std::cout << "m1=" << m << "\n"
		<< "m2=" << polyset_type::toMatrix( polyset_type::toType( m ) ) << "\n"
		<< ublas::norm_frobenius( polyset_type::toMatrix( polyset_type::toType( m ) ) - m ) << "\n";
    LIFE_ASSERT( ublas::norm_frobenius( polyset_type::toMatrix( polyset_type::toType( m ) ) -
					m ) < 1e-10 )( m ).warn ( "invalid transformation" );
    this->setCoefficient( polyset_type::toType( m ), true );
  }

  OrthonormalPolynomialSet<Dim, Order, Scalar,T, Simplex > toScalar() const
  {
    return OrthonormalPolynomialSet<Dim, Order, Scalar,T, Simplex >();
  }

  /**
   * \return the family name of the polynomial set
   */
  std::string familyName() const { return "dubiner"; }


  points_type points() const { return points_type(); }
  points_type points( int f ) const { return points_type(); }
};

  template<uint16_type Dim,
	   uint16_type Order,
	   template<uint16_type> class PolySetType,
	   typename T>
  const uint16_type OrthonormalPolynomialSet<Dim, Order,PolySetType,T, Simplex>::nLocalDof;


template<uint16_type Dim,
         uint16_type Order,
         template<uint16_type> class PolySetType,
         typename T>
class OrthonormalPolynomialSet<Dim, Order, PolySetType, T, SimplexProduct>
  :
    public PolynomialSet<Legendre<Dim, Order, Normalized<true>, T>, PolySetType >
{
  typedef PolynomialSet<Legendre<Dim, Order, Normalized<true>, T>, PolySetType > super;
public:

  static const uint16_type nDim = Dim;
  static const uint16_type nOrder = Order;

  typedef OrthonormalPolynomialSet<Dim, Order, PolySetType, T, SimplexProduct> self_type;
  typedef self_type component_basis_type;

  typedef typename super::polyset_type polyset_type;
  static const bool is_tensor2 = polyset_type::is_tensor2;
  static const bool is_vectorial = polyset_type::is_vectorial;
  static const bool is_scalar = polyset_type::is_scalar;
  static const bool is_continuous = false;
  static const bool is_modal = true;
  static const uint16_type nComponents = polyset_type::nComponents;
  typedef typename super::component_type component_type;
  typedef T value_type;
  typedef Legendre<Dim, Order, Normalized<true>, T> basis_type;
  typedef SimplexProduct<Dim, Order, Dim> convex_type;
  typedef typename matrix_node<value_type>::type points_type;

  template<int O>
  struct convex
  {
    typedef SimplexProduct<Dim, O, Dim> type;
  };
  typedef Reference<convex_type, nDim, nOrder, nDim, value_type> reference_convex_type;

  typedef typename super::polynomial_type polynomial_type;

  //!< Number of degrees of freedom per vertex
  static const uint16_type nDofPerVertex = convex_type::nbPtsPerVertex;
  //!< Number of degrees  of freedom per edge
  static const uint16_type nDofPerEdge = convex_type::nbPtsPerEdge;
  //!< Number of degrees  of freedom per face
  static const uint16_type nDofPerFace = convex_type::nbPtsPerFace;

  //!< Number of degrees  of freedom per volume
  static const uint16_type nDofPerVolume = convex_type::nbPtsPerVolume;

  static const uint16_type nLocalDof = convex_type::numPoints;

  static const uint16_type nDof = nLocalDof;
  static const uint16_type nNodes = nDof;
  static const uint16_type nDofGrad = super::nDim*nDof;
  static const uint16_type nDofHess = super::nDim*super::nDim*nDof;

  OrthonormalPolynomialSet()
    :
    super( basis_type() )

  {
    ublas::matrix<value_type> m( ublas::identity_matrix<value_type>( nComponents*convex_type::polyDims( nOrder ) ) );
    if ( is_tensor2 )
      std::cout << "[orthonormalpolynomialset] m = " << m << "\n";
    LIFE_ASSERT( ublas::norm_frobenius( polyset_type::toMatrix( polyset_type::toType( m ) ) -
					m ) < 1e-10 )( m ).warn ( "invalid transformation" );
    this->setCoefficient( polyset_type::toType( m ), true );
  }

  OrthonormalPolynomialSet<Dim, Order, Scalar,T, SimplexProduct > toScalar() const
  {
    return OrthonormalPolynomialSet<Dim, Order, Scalar,T, SimplexProduct >();
  }
  std::string familyName() const { return "legendre"; }
  points_type points() const { return points_type(); }
  points_type points( int f ) const { return points_type(); }
};

  template<uint16_type Dim,
	   uint16_type Order,
	   template<uint16_type> class PolySetType,
	   typename T>
  const uint16_type OrthonormalPolynomialSet<Dim, Order,PolySetType,T, SimplexProduct>::nLocalDof;



} // Life
#endif /* __PolynomialSet_H */

