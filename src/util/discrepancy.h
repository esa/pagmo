/*****************************************************************************
 *   Copyright (C) 2004-2013 The PaGMO development team,                     *
 *   Advanced Concepts Team (ACT), European Space Agency (ESA)               *
 *   http://apps.sourceforge.net/mediawiki/pagmo                             *
 *   http://apps.sourceforge.net/mediawiki/pagmo/index.php?title=Developers  *
 *   http://apps.sourceforge.net/mediawiki/pagmo/index.php?title=Credits     *
 *   act@esa.int                                                             *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program; if not, write to the                           *
 *   Free Software Foundation, Inc.,                                         *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.               *
 *****************************************************************************/

#ifndef PAGMO_UTIL_DISCREPANCY_H
#define PAGMO_UTIL_DISCREPANCY_H

#include <iostream>
#include <vector>
#include <math.h>
#include <algorithm>
#include <boost/shared_ptr.hpp>

#include "../config.h"
#include "../exceptions.h"

namespace pagmo{ namespace util {

/// Discrepancy namespace.
/**
 * Utilities to generate low-discrepancy sequences in the
 * hypercube or in a simplex
 *
 * @see http://en.wikipedia.org/wiki/Quasi-Monte_Carlo_method
*/
namespace discrepancy {

//! @cond
double van_der_corput(unsigned int n, unsigned int base);
unsigned int prime ( int n );
unsigned int prime_ge ( unsigned int n );
//! @endcond

class base;
typedef boost::shared_ptr<base> base_ptr;

class __PAGMO_VISIBLE base
{
public:
	base(unsigned int dim, unsigned int count = 1) : m_dim(dim), m_count(count) {}
	virtual std::vector<double> operator()() = 0;
	virtual std::vector<double> operator()(size_t) = 0;
	virtual base_ptr clone() const = 0;
	virtual~base() {}
protected:
	unsigned int m_dim;
	unsigned int m_count;
};

/// Halton quasi-random point sequence
/**
 * Class that generates a quasi-random sequence of
 * points in the unit hyper cube using the Halton sequence
 *
 * @see http://en.wikipedia.org/wiki/Halton_sequence
 * @author dario.izzo@gmail.com
*/
class __PAGMO_VISIBLE halton : public base
{
	public:
		/// Constructor
		/**
		 * @param[in] dim dimension of the hypercube
		 * @param[in] count starting point of the sequence (the first point is  [0.5,0.33333, ....])
		 *
		 * @throws value_error if dim>10
		*/
		halton(unsigned int dim, unsigned int count = 1) : base(dim,count), m_primes() {
			if (dim >10) {
				pagmo_throw(value_error,"Halton sequences should not be used in dimension >10");
			}
			for (size_t i=1; i<=dim; ++i) {
				m_primes.push_back(prime(i));
			}
		}
		/// Clone method.
		base_ptr clone() const
		{
			return base_ptr(new halton(*this));
		}
		/// Operator ()
		/**
		 * Returns the next point in the sequence
		 *
		 * @return an std::vector<double> containing the next point
		 */
		std::vector<double> operator()() {
			std::vector<double> retval;
			for (size_t i=0; i<m_dim; ++i) {
				retval.push_back(van_der_corput(m_count,m_primes.at(i)));
			}
			m_count++;
			return retval;
		}
		/// Operator (size_t n)
		/**
		 * Returns the n-th point in the sequence
		 *
		 * @param[in] n the point along the sequence to be returned
		 * @return an std::vector<double> containing the n-th point
		 */
		std::vector<double> operator()(size_t n) {
			if (n == 0) {
				pagmo_throw(value_error,"Halton sequence first point id is 1");
			}
			std::vector<double> retval;
			for (size_t i=0; i<m_primes.size(); ++i) {
				retval.push_back(van_der_corput(n,m_primes.at(i)));
			}
			m_count = n+1;
			return retval;
		}
	private:
		std::vector<unsigned int> m_primes;
};


/// Faure quasi-random point sequence
/**
 * Class that generates a quasi-random sequence of
 * points in the unit hyper cube using the Faure sequence.
 * The code wraps original routines from the link below.
 *
 * @see http://people.sc.fsu.edu/~jburkardt/cpp_src/faure/faure.html
 * @author dario.izzo@gmail.com
*/
class __PAGMO_VISIBLE faure : public base
{
	public:
	/// Constructor
	/**
	 * @param[in] dim dimension of the hypercube
	 * @param[in] count starting point of the sequence
	 *
	 * @throws value_error if dim not in [2,23]
	*/
	faure(unsigned int dim, unsigned int count = 1) : base(dim, count), m_coef(NULL), m_hisum_save(-1), m_qs(-1), m_ytemp(NULL) {
			if (dim >23 || dim <2) {
				pagmo_throw(value_error,"Faure sequences can have dimension [2,23]");
			}
		}
	/// Clone method.
	base_ptr clone() const
	{
		return base_ptr(new faure(*this));
	}
	/// Operator ()
	/**
	 * Returns the next point in the sequence
	 *
	 * @return an std::vector<double> containing the next point
	 */
	std::vector<double> operator()() {
		std::vector<double> retval(m_dim,0.0);
		faure_orig(m_dim, &m_count, &retval[0]);
		return retval;
	}
	/// Operator (size_t n)
	/**
	 * Returns the n-th point in the sequence
	 *
	 * @param[in] n the point along the sequence to be returned
	 * @return an std::vector<double> containing the n-th point
	 */
	std::vector<double> operator()(size_t n) {
		m_count = n;
		std::vector<double> retval(m_dim,0.0);
		faure_orig(m_dim, &m_count, &retval[0]);
		return retval;
	}
	private:
		int *binomial_table ( int qs, int m, int n );
		void faure_orig ( unsigned int dim_num, unsigned int *seed, double quasi[] );
		double *faure_generate ( int dim_num, int n, int skip );
		int i4_log_i4 ( int i4, int j4 );
		int i4_min ( int i1, int i2 );
		int i4_power ( int i, int j );
	private:
		int *m_coef;
		int m_hisum_save;
		int m_qs;
		int *m_ytemp;

};

class __PAGMO_VISIBLE project_2_simplex
{
	public:
		project_2_simplex(unsigned int dim) : m_dim(dim) {}
		std::vector<double> operator()(std::vector<double> retval) const {
			if (retval.size() == m_dim-1) {
				std::sort(retval.begin(),retval.end());
				retval.insert(retval.begin(),0.0);
				retval.push_back(1.0);
				double cumsum=0;
				for (unsigned int i = 0; i<retval.size()-1;++i) {
					retval[i] = retval[i+1] - retval[i];
					cumsum += retval[i];
				}
				retval.pop_back();
				for (unsigned int i = 0; i<retval.size();++i) {
					retval[i] /= cumsum;
				}
				return retval;
			}
			else {
				pagmo_throw(value_error,"To project on this simplex you need a point in dimension m_dim-1");
			}
		}
	private:
		unsigned int m_dim;
};

class __PAGMO_VISIBLE simplex : public base
{
public:
	simplex(unsigned int dim, unsigned int count) : base(dim,count), m_generator(dim-1), m_projector(dim) {}
	/// Clone method.
	base_ptr clone() const
	{
		return base_ptr(new simplex(*this));
	}
	std::vector<double> operator()() {
		std::vector<double> tmp = m_generator();
//std::cout << tmp[0] << " " << tmp[1] << " " << std::endl;
		std::vector<double> retval = m_projector(tmp);
		return retval;
	}
	std::vector<double> operator()(size_t n) {
		std::vector<double> retval = m_projector(m_generator(n));
		return retval;
	}
private:
	halton m_generator;
	project_2_simplex m_projector;
};

}}} //namespace discrepancy

#endif
