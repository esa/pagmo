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
 * hypercube or on a simplex
 *
 * @see http://en.wikipedia.org/wiki/Quasi-Monte_Carlo_method
*/
namespace discrepancy {

//! @cond
double van_der_corput(unsigned int n, unsigned int base);
unsigned int prime ( int n );
unsigned int prime_ge ( unsigned int n );
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
//! @endcond

class base;
typedef boost::shared_ptr<base> base_ptr;

//---------------------------------------------------------
//--------------------BASE CLASS---------------------------
//---------------------------------------------------------

/// Base low-discrepancy sequence class
/**
 * This class cannot be instantiated as it contains pure virtual members. All classes
 * that generate quasi-random sequences with low-discrepancy must inherit from this class. 
 * 
 * @author Dario Izzo (dario.izzo@gmail.com)
 */
class __PAGMO_VISIBLE base
{
public:
	/// Constructor
	/**
	 * @param[in] dim hypercube dimension
	 * @param[in] count starting point of the sequence
	*/
	base(unsigned int dim, unsigned int count = 1) : m_dim(dim), m_count(count) {}
	/// Operator ()
	/**
	 * Returns the next point in the sequence. Must be implemented in the derived class
	 *
	 * @return an std::vector<double> containing the next point
	 */
	virtual std::vector<double> operator()() = 0;
	/// Operator (size_t n)
	/**
	 * Returns the n-th point in the sequence. Must be implemented in the derived class
	 *
	 * @param[in] n the point along the sequence to be returned
	 * @return an std::vector<double> containing the n-th point
	 */
	virtual std::vector<double> operator()(unsigned int n) = 0;
	/// Clone method for dynamic polymorphism
	virtual base_ptr clone() const = 0;
	/// Virtual destructor. Required as the class contains pure virtual methods
	virtual~base();
protected:
	unsigned int m_dim;
	unsigned int m_count;
};

//---------------------------------------------------------
//--------------------DERIVED CLASSES----------------------
//---------------------------------------------------------

/// Halton quasi-random point sequence
/**
 * Class that generates a quasi-random sequence of
 * points in the unit hyper-cube using the Halton sequence
 *
 * @see http://en.wikipedia.org/wiki/Halton_sequence
 * @author dario.izzo@gmail.com
*/
class __PAGMO_VISIBLE halton : public base
{
	public:
		halton(unsigned int dim, unsigned int count = 1);
		base_ptr clone() const;
		std::vector<double> operator()();
		std::vector<double> operator()(unsigned int n);
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
	faure(unsigned int dim, unsigned int count = 1);
	base_ptr clone() const;
	std::vector<double> operator()();
	std::vector<double> operator()(unsigned int n);
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

/// Halton sequence projected on a simplex
/**
 * Class that generates a quasi-random sequence of
 * points on a n-dimensional simplex. In essence we sample a point such that:
 * \f$ \sum_i x_{i} = 1 \f$
 * And we do this using Halton sequence. The algorithm is original (as far as we know).
 *
 * @author dario.izzo@gmail.com
*/
class __PAGMO_VISIBLE simplex : public base
{
public:
	simplex(unsigned int dim, unsigned int count);
	base_ptr clone() const;
	std::vector<double> operator()();
	std::vector<double> operator()(unsigned int n);
private:
	halton m_generator;
	project_2_simplex m_projector;
};

}}} //namespace discrepancy

#endif
