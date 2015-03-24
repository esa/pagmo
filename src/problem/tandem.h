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

#ifndef PAGMO_PROBLEM_TANDEM_H
#define PAGMO_PROBLEM_TANDEM_H

#include <string>
#include <vector>

#include "../config.h"
#include "../serialization.h"
#include "../types.h"
#include "../AstroToolbox/mga_dsm.h"
#include "../AstroToolbox/misc4Tandem.h"
#include "base.h"

namespace pagmo{ namespace problem {

/// TandEM problem
/**
 * \image html tandem.jpg "Encedalus orbiting around Saturn"
 * \image latex tandem.jpg "Encedalus orbiting around Saturn" width=5cm
 *
 * TandEM primary goals are to understand the atmosphere, surface and interior,
 * to determine the chemistry, and to derive constraints on the origin and evolution of
 * Titan and of the saturnian system as a whole, with an emphasis on Enceladus.
 *
 * TandEM had been proposed as an L-class candidate mission for the Cosmic Vision 2015-2025
 * programme in response to the Call for Proposals issued by ESA in March 2007,
 * and has been one of the mission concepts selected by the SSAC in October 2007
 * for an Assessment study. The SSAC had decided that an Outer Planet mission should be one
 * of the L-class mission candidates for the first slice of the Cosmic Vision plan,
 * and recommended that both TandEM and Laplace (a mission to the Jupiter system)
 * concepts should be studied initially, with a decision on which one to pursue to take place
 * following a better definition of the mission's characteristics.

 * In early February 2009 ESA and NASA jointly announced that the Europa-Jupiter System Mission,
 * or Laplace (see problem::laplace), would be the candidate for the first L mission.
 *
 * We transcribe here TandEM interplanetary trajectory transfer problem as an MGA-DSM
 * problem and we allow to instanciate 25 different types, depending on the fly-by sequence adopted.
 * The mission data are taken from those used by a joint ESA/NASA working group that performed, together with industrial partners
 * a preliminary trajectory design for the TandEM mission during the first months of 2008.
 * Please refer to http://www.esa.int/gsp/ACT/inf/op/globopt/TandEM.htm to select the proper instance. A default
 * choice is 6, which corrsponds to an EVEES fly-by sequence. The possibility of adding one additional constraint
 * on the time-of-flight is given. If such a constraint is specified the components x[4]-x[7] are no longer time of flights
 * but they represent the percentage of the remaining time of flight to be used in one leg. Thus the linear constraint
 * on the time of flight is automatically transformed into a box constraint.
 *
 * The problem is transcribed as an MGA-DSM problem allowing one chemical manouvre per trajectory leg.
 * The objective function is defined as \f$ f = -\log(m_f)\f$ where \f$m_f\f$ is the final spacecraft mass.
 * The problem is also part of the Global Trajectory Optimization database (GTOP)
 *
 * tandem is a box constrained single objective, continuous optimization problem of dimension 18.
 *
 * @see http://www.esa.int/gsp/ACT/inf/op/globopt/TandEM.htm
 * @author Dario Izzo (dario.izzo@esa.int)
 */
class __PAGMO_VISIBLE tandem: public base
{
	public:
		tandem(const int problemid = 6, const double tof_ = -1);
		base_ptr clone() const;
		std::string pretty(const std::vector<double> &x) const;
		std::string get_name() const;
	protected:
		void objfun_impl(fitness_vector &, const decision_vector &) const;
		void set_sparsity(int &, std::vector<int> &, std::vector<int> &) const;
	private:
		friend class boost::serialization::access;
		template <class Archive>
		void serialize(Archive &ar, const unsigned int)
		{
			ar & boost::serialization::base_object<base>(*this);
			ar & problem;
			ar & const_cast<double &>(tof);
			ar & copy_of_x;
		}
		static const int Data[24][5]; // [DS] These two arrays are not serialized as they 
		static const int sequence[5]; // are declared as static consts
		mgadsmproblem problem;
		const double tof;
		mutable std::vector<double> copy_of_x;

};

}}

BOOST_CLASS_EXPORT_KEY(pagmo::problem::tandem)

#endif // PAGMO_PROBLEM_TANDEM_H
