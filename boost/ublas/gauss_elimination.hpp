/**
	Boost.uBLAS Solver algorithms competency test for GSoC15
	Author : Ganesh Prasad Sahoo
	March 2015
*/

#ifndef BOOST_UBLAS_GAUSS_ELIMINATION_HPP_
#define BOOST_UBLAS_GAUSS_ELIMINATION_HPP_

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/lu.hpp>

#include <algorithm>

// Gaussian Elimination algorithm

namespace boost{ namespace numeric{ namespace ublas{
	
	// gaussian elimination without pivoting
	// remarks:
	// fast but numerically unstable
	template<typename M>
	typename M::size_type gaussian_elimination_no_pivot(M &m){
		typedef M matrix_type;
		typedef typename M::size_type size_type;
		typedef typename M::value_type value_type;

		size_type singular = 0;
		size_type size1 = m.size1();
		size_type size2 = m.size2();
		size_type size = (std::min) (size1, size2);

#if BOOST_UBLAS_TYPE_CHECK
		matrix_type cm(m);
		matrix_type l(size, size); // matrix_type l(identity_matrix <value_type> (size)); does not work here ??? why ??
		l(size-1, size-1) = 1;
#endif

		for(size_type k=0; k<size-1; k++){
#if BOOST_UBLAS_TYPE_CHECK
			l(k, k) = 1;
#endif
			if(m(k, k) != value_type(0)){		
				for(size_type j=k+1; j<size; j++){
					value_type coeff = m(j, k) / m(k, k);
#if BOOST_UBLAS_TYPE_CHECK
					l(j, k) = coeff;
#endif
					matrix_row<M> m_j (row(m, j));
					matrix_row<M> m_k (row(m, k));
					project(m_j, range(k, size)) -= coeff * project(m_k, range(k, size));
				}
			}
			else if( singular == 0 ){
				singular = k+1;
			}
		}

#if BOOST_UBLAS_TYPE_CHECK
		BOOST_UBLAS_CHECK ( singular != 0 ||
							detail::expression_type_check ( prod( triangular_adaptor<matrix_type, unit_lower> (l),
																 triangular_adaptor<matrix_type, upper> (m) ),
															cm), internal_logic());
#endif
		return singular;
	}



	
	namespace gauss_aux{

		// this is practically useless, boost::range::max_element should be used, 
		// but my implementation went wrong somewhere
		// hence this retreat to std::max_element

		template<typename T>
		BOOST_UBLAS_INLINE
		bool compare_abs(const T &a, const T &b){
			return std::abs(a) > std::abs(b);
		}

		template<typename MC>
		BOOST_UBLAS_INLINE
		typename MC::size_type gauss_pivot_argmax(MC &mc, typename MC::size_type k){
			typedef MC column_type;
			typedef typename MC::size_type size_type;
			typedef typename MC::value_type value_type;
			typedef typename MC::iterator iterator_type;
		
			iterator_type pos = std::max_element(mc.begin()+k, mc.end(), compare_abs<value_type>);
			size_type result = pos - mc.begin();

			return result;
		}
	}

	// gauss elimination with partial pivoting
	// without permutation matrix
	// remarks:
	// relatively slow but numerically stable

	template<typename M>
	typename M::size_type gaussian_elimination(M &m){
		typedef M matrix_type;
		typedef typename M::size_type size_type;
		typedef typename M::value_type value_type;
		typedef vector<value_type> vector_type;
		
		size_type singular = 0;
		size_type size1 = m.size1();
		size_type size2 = m.size2();
		size_type size = std::min (size1, size2);

#if BOOST_UBLAS_TYPE_CHECK
		matrix_type cm(m);
		identity_matrix<value_type> dummy(size);
		matrix_type l(dummy); // matrix_type l(identity_matrix<value_type>(size)); does not work
#endif

		for(size_type k=0; k<size-1; k++){
			matrix_column<matrix_type> m_k(m, k);
			size_type i_max = gauss_aux::gauss_pivot_argmax(m_k, k);
			if( m(i_max, k) != value_type(0)){
				row(m, k).swap(row(m, i_max));
#if BOOST_UBLAS_TYPE_CHECK
				row(l, k).swap(row(l, i_max));
#endif
				for(size_type i = k+1; i < size; i++){
					BOOST_UBLAS_CHECK ( m(k, k) != 0, divide_by_zero() );
					value_type coeff = m(i, k) / m(k, k);
#if BOOST_UBLAS_TYPE_CHECK
					l(i, k) = coeff;
#endif
					matrix_row<M> m_i (row(m, i));
					matrix_row<M> m_k (row(m, k));
					project(m_i, range(k, size)) -= coeff * project(m_k, range(k, size));

					m(i, k) = 0;
				}
			}
			else if(singular == 0){
				singular = k+1;
			}
		}
#if BOOST_UBLAS_TYPE_CHECK
		BOOST_UBLAS_CHECK ( singular != 0 ||
							detail::expression_type_check ( prod( triangular_adaptor<matrix_type, unit_lower> (l),
																 triangular_adaptor<matrix_type, upper> (m) ),
															cm), internal_logic());
#endif
		return singular;
	}






}}}

#endif
