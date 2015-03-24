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

#ifndef PAGMO_PROBLEM_CEC2013_H
#define PAGMO_PROBLEM_CEC2013_H

#include <string>

#include "../serialization.h"
#include "../types.h"
#include "base.h"

namespace pagmo{ namespace problem {

/// The CEC 2013 problems: Real-Parameter Single Objective Optimization Competition
/**
 *
 * This class allows to instantiate any of the 28 problems of the competition on real-parameter
 * single objective optimization problems that was organizd in the framework of the
 * 2013 IEEE Congress on Evolutionary Computation.
 *
 * NOTE: It requires the data files that can be downloaded from the link below
 * upon construction, it expects to find in the folder indicated by
 * the constructor argument std::string for two files named M_Dxx.txt and shift_data.txt.
 *
 * NOTE 2: all problems are unconstrained continuous single objective problems.
 *
 * @see http://www.ntu.edu.sg/home/EPNSugan/index_files/CEC2013/CEC2013.htm
 *
 * @author Dario Izzo (dario.izzo@gmail.com)
 */

class __PAGMO_VISIBLE cec2013 : public base
{
	public:
	cec2013(unsigned int = 1, problem::base::size_type = 30, const std::string & = "input_data/");
		base_ptr clone() const;
		std::string get_name() const;

		/** @name Getters.*/
		//@{
		/// Returns the origin shift used by the problem
		/**
		 * @returns the origin shift
		 *
		 */
		std::vector<double> origin_shift() const {return m_origin_shift;}
		//@}
	protected:
		void objfun_impl(fitness_vector &, const decision_vector &) const;
	private:
		void sphere_func (const double *, double *, int , const double *,const double *, int) const; /* Sphere */
		void ellips_func(const double *, double *, int , const double *,const double *, int) const; /* Ellipsoidal */
		void bent_cigar_func(const double *, double *, int , const double *,const double *, int) const; /* Discus */
		void discus_func(const double *, double *, int , const double *,const double *, int) const;  /* Bent_Cigar */
		void dif_powers_func(const double *, double *, int , const double *,const double *, int) const;  /* Different Powers */
		void rosenbrock_func (const double *, double *, int , const double *,const double *, int) const; /* Rosenbrock's */
		void schaffer_F7_func (const double *, double *, int , const double *,const double *, int) const; /* Schwefel's F7 */
		void ackley_func (const double *, double *, int , const double *,const double *, int) const; /* Ackley's */
		void rastrigin_func (const double *, double *, int , const double *,const double *, int) const; /* Rastrigin's  */
		void weierstrass_func (const double *, double *, int , const double *,const double *, int) const; /* Weierstrass's  */
		void griewank_func (const double *, double *, int , const double *,const double *, int) const; /* Griewank's  */
		void schwefel_func (const double *, double *, int , const double *,const double *, int) const; /* Schwefel's */
		void katsuura_func (const double *, double *, int , const double *,const double *, int) const; /* Katsuura */
		void bi_rastrigin_func (const double *, double *, int , const double *,const double *, int) const; /* Lunacek Bi_rastrigin */
		void grie_rosen_func (const double *, double *, int , const double *,const double *, int) const; /* Griewank-Rosenbrock  */
		void escaffer6_func (const double *, double *, int , const double *,const double *, int) const; /* Expanded Scaffer¡¯s F6  */
		void step_rastrigin_func (const double *, double *, int , const double *,const double *, int) const; /* Noncontinuous Rastrigin's  */
		void cf01 (const double *, double *, int , const double *,const double *, int) const; /* Composition Function 1 */
		void cf02 (const double *, double *, int , const double *,const double *, int) const; /* Composition Function 2 */
		void cf03 (const double *, double *, int , const double *,const double *, int) const; /* Composition Function 3 */
		void cf04 (const double *, double *, int , const double *,const double *, int) const; /* Composition Function 4 */
		void cf05 (const double *, double *, int , const double *,const double *, int) const; /* Composition Function 5 */
		void cf06 (const double *, double *, int , const double *,const double *, int) const; /* Composition Function 6 */
		void cf07 (const double *, double *, int , const double *,const double *, int) const; /* Composition Function 7 */
		void cf08 (const double *, double *, int , const double *,const double *, int) const; /* Composition Function 8 */
		void shiftfunc (const double*,double*,int,const double*) const;
		void rotatefunc (const double *,double*,int, const double *) const;
		void asyfunc (const double *, double *, int, double) const;
		void oszfunc (const double *, double *, int) const;
		void cf_cal(const double *, double *, int, const double *,double *,double *,double *,int) const;

		friend class boost::serialization::access;
		template <class Archive>
		void serialize(Archive &ar, const unsigned int)
		{
			ar & boost::serialization::base_object<base>(*this);
			ar & const_cast<unsigned int&>(m_problem_number);
			ar & m_rotation_matrix;
			ar & m_origin_shift;
		}
	const unsigned int m_problem_number;
	std::vector<double> m_rotation_matrix;
	std::vector<double> m_origin_shift;

	// These are pre-allocated for speed, need not to be serialized
	mutable std::vector<double> m_y;
	mutable std::vector<double> m_z;
};

}} //namespaces

BOOST_CLASS_EXPORT_KEY(pagmo::problem::cec2013)

#endif
