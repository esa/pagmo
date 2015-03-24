/*****************************************************************************
 *   Copyright (C) 2004-2015 The PaGMO development team,                     *
 *   Advanced Concepts Team (ACT), European Space Agency (ESA)               *
 *                                                                           *
 *   https://github.com/esa/pagmo                                            *
 *                                                                           *
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
#include "../rng.h"

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
		std::vector<double> operator()(std::vector<double> retval) const;
	private:
		unsigned int m_dim;
};
//! @endcond

class base;
/// Smart pointer to the base discrepancy class
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
	/// Hypercube dimension where sampling with low-discrepancy
	unsigned int m_dim;
	/// Starting point of the sequence (can be used to skip initial values)
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


/// Sobol quasi-random point sequence
/**
 * Class that generates a quasi-random sequence of
 * points in the unit hyper cube using the Sobol sequence.
 * The code wraps original routines from the link below.
 *
 * @see http://people.sc.fsu.edu/~jburkardt/cpp_src/sobol/sobol.html
 * @author c.ortega.absil@gmail.com
*/
class __PAGMO_VISIBLE sobol : public base
{
	public:
		sobol(unsigned int dim, unsigned int count);
		base_ptr clone() const;
		std::vector<double> operator()();
		std::vector<double> operator()(unsigned int n);
	private:
		int i8_bit_lo0 ( long long int n );
		void i8_sobol ( unsigned int dim_num, long long int *seed, double quasi[ ] );
	private:
		unsigned int m_dim_num_save;
		bool m_initialized;
		long long int m_maxcol;
		long long int m_seed_save;
		double recipd;
		long long int lastq[1111]; //1111 is maximum dimension.
		long long int poly[1111];
		long long int v[1111][62]; //2^62 is approx. limit of points requested.
};

/// Latin Hypercube Sampling
/**
 * Class that generates a latin hypersquare sampling
 * in the unit hyper cube.
 * The code wraps original routines from the link below.
 *
 * @see http://people.sc.fsu.edu/~jburkardt/cpp_src/latin_random/latin_random.html
 * @author c.ortega.absil@gmail.com
*/


class __PAGMO_VISIBLE lhs : public base
{
	public:
		lhs(unsigned int dim, unsigned int count);
		base_ptr clone() const;
		std::vector<double> operator()();
		std::vector<double> operator()(unsigned int n);
	private:
		std::vector<double> latin_random ( unsigned int dim_num, unsigned int point_num);
		unsigned int *perm_uniform ( unsigned int n);
	private:
		bool m_initialised;
		std::vector<double> m_set;
		unsigned int m_next;
};

}}} //namespace discrepancy

#endif
