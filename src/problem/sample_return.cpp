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

#include <string>

#include "sample_return.h"
#include <keplerian_toolbox/core_functions/convert_dates.h>



namespace pagmo { namespace problem {

static const int EA[2] = {3,10};
static const int AE[2] = {10,3};

/// Constructor
/**
 * Instantiates the sample_return problem
 */
sample_return::sample_return(const ::kep_toolbox::planet::base &asteroid, const double &Tmax):base(12), m_target(asteroid.clone()),
	 m_leg1(total_DV_rndv,EA,2,0,0,0,0,0), m_leg2(total_DV_rndv,AE,2,0,0,0,0,0),x_leg1(6),x_leg2(6),m_Tmax(Tmax)
{
	::kep_toolbox::epoch start_lw(2020,1,1);
	::kep_toolbox::epoch end_lw(2035,1,1);
	const double lb[12] = {start_lw.mjd2000(),0  ,0,0,0.1  ,0.001, 20 ,0,0,0,0.1   ,0.001};
	const double ub[12] = {end_lw.mjd2000()  ,5.5,1,1,0.7, 0.999, 50  ,6,1,1,1     ,0.999};
	set_bounds(lb,lb+12,ub,ub+12);

	for (int i = 0;i<6;++i)
	{
		m_leg1.asteroid.keplerian[i] = m_target->compute_elements(::kep_toolbox::epoch(0))[i];
		m_leg2.asteroid.keplerian[i] = m_target->compute_elements(::kep_toolbox::epoch(0))[i];
	}

	//convert to JPL unit format .... (check mgadsm)
	m_leg1.asteroid.keplerian[0] /= ASTRO_AU;
	m_leg2.asteroid.keplerian[0] /= ASTRO_AU;
	for (int i = 2;i<6;++i)
	{
		m_leg1.asteroid.keplerian[i] *= ASTRO_RAD2DEG;
		m_leg2.asteroid.keplerian[i] *= ASTRO_RAD2DEG;
	}

	m_leg1.asteroid.epoch = ::kep_toolbox::mjd20002mjd(0);
	m_leg2.asteroid.epoch = ::kep_toolbox::mjd20002mjd(0);
}

/// Clone method.
base_ptr sample_return::clone() const
{
	return base_ptr(new sample_return(*this));
}

/// Implementation of the objective function.
void sample_return::objfun_impl(fitness_vector &f, const decision_vector &x) const
{
	//We split the decision vector in the two legs
	std::copy(x.begin(),x.begin()+6,x_leg1.begin());
	std::copy(x.begin()+6,x.begin()+12,x_leg2.begin());

	x_leg1[4] = x_leg1[4] * m_Tmax;
	x_leg2[4] = (m_Tmax - x_leg1[4] - x_leg2[0]) * x_leg2[4];

	//We account for the waiting time
	x_leg2[0] += x_leg1[0] + x_leg1[4];
	double dummy = 0;
	MGA_DSM(x_leg1, m_leg1, dummy);
	MGA_DSM(x_leg2, m_leg2, dummy);
	f[0] = m_leg1.DV[0] + m_leg1.DV[1] + m_leg1.DV[2] +
		m_leg2.DV[0] + m_leg2.DV[1] + std::max(0.0,m_leg2.DV[2] - 5.5);

}

/// Outputs a stream with the trajectory data
/**
 * While the chromosome contains all necessary information to describe a trajectory, mission analysits
 * often require a different set of data to evaluate its use. This method outputs a stream with
 * information on the trajectory that is otherwise 'hidden' in the chromosome
 *
 * \param[in] x chromosome representing the trajectory in the optimization process
 * \returns an std::string with launch dates, DV magnitues and other information on the trajectory
 */

std::string sample_return::pretty(const std::vector<double> &x) const
{
	std::ostringstream s;
	s.precision(15);
	s << std::scientific;
	//We split the decision vector in the two legs
	std::copy(x.begin(),x.begin()+6,x_leg1.begin());
	std::copy(x.begin()+6,x.begin()+12,x_leg2.begin());

	x_leg1[4] = x_leg1[4] * m_Tmax;
	x_leg2[4] = (m_Tmax - x_leg1[4] - x_leg2[0]) * x_leg2[4];

	//We account for the waiting time
	x_leg2[0] += x_leg1[0] + x_leg1[4];
	double dummy = 0;
	MGA_DSM(x_leg1, m_leg1, dummy);
	MGA_DSM(x_leg2, m_leg2, dummy);

	s << "Departure epoch (mjd2000):\t" << x[0] << '\n';
	s << "Departure epoch:\t\t" << ::kep_toolbox::epoch(x[0],::kep_toolbox::epoch::MJD2000) << '\n';
	s << "Escape velocity:\t\t" << m_leg1.DV[0] << '\n';
	s << "dsm1 epoch:\t\t\t" << ::kep_toolbox::epoch(x[0] + x[5]*x[4],::kep_toolbox::epoch::MJD2000) << '\n';
	s << "dsm1 magnitude\t\t\t" << m_leg1.DV[1] << " \n";
	s << "Asteroid arrival epoch: \t" << ::kep_toolbox::epoch(x[0] + x[4],::kep_toolbox::epoch::MJD2000) << '\n';
	s << "Breaking manouvre:\t\t" << m_leg1.DV[2] << '\n' << std::endl;

	s << "Departure epoch (mjd2000):\t" << x[0] + x[4] + x[6] << '\n';
	s << "Departure epoch:\t\t" << ::kep_toolbox::epoch(x[0] + x[4] + x[6],::kep_toolbox::epoch::MJD2000) << '\n';
	s << "Escape velocity:\t\t" << m_leg2.DV[0] << '\n';
	s << "dsm2 epoch:\t\t\t" << ::kep_toolbox::epoch(x[0]+x[4]+x[6]+x[11]*x[10],::kep_toolbox::epoch::MJD2000) << '\n';
	s << "dsm2 magnitude\t\t\t" << m_leg2.DV[1] << " \n";
	s << "Earth arrival epoch: \t\t" << ::kep_toolbox::epoch(x[0]+x[4]+x[6]+x[10],::kep_toolbox::epoch::MJD2000) << '\n';
	s << "Arrival Vinf:\t\t\t" << m_leg2.DV[2] << '\n';
	s << "Total time of flight:\t\t" << x[4]+x[6]+x[10] << '\n' << std::endl;

	s << "Earth-ephemerides at departure:\t\t" << m_leg1.r[0][0] << " " << m_leg1.r[0][1] << " " << m_leg1.r[0][2] << std::endl;
	s << "Asteroid-ephemerides at arrival:\t" << m_leg1.r[1][0] << " " << m_leg1.r[1][1] << " " << m_leg1.r[1][2] << std::endl;
	s << "Asteroid-ephemerides at departure:\t" << m_leg2.r[0][0] << " " << m_leg2.r[0][1] << " " << m_leg2.r[0][2] << std::endl;
	s << "Earth-ephemerides at arrival:\t\t" << m_leg2.r[1][0] << " " << m_leg2.r[1][1] << " " << m_leg2.r[1][2] << std::endl;



	return s.str();
}

/// Computes the Delta-Vs
std::vector<double> sample_return::get_delta_v(const std::vector<double> &x) const {
	//We split the decision vector in the two legs
	std::copy(x.begin(),x.begin()+6,x_leg1.begin());
	std::copy(x.begin()+6,x.begin()+12,x_leg2.begin());

	x_leg1[4] = x_leg1[4] * m_Tmax;
	x_leg2[4] = (m_Tmax - x_leg1[4] - x_leg2[0]) * x_leg2[4];

	//We account for the waiting time
	x_leg2[0] += x_leg1[0] + x_leg1[4];
	double dummy = 0;
	MGA_DSM(x_leg1, m_leg1, dummy);
	MGA_DSM(x_leg2, m_leg2, dummy);
	std::vector<double> retval;
	for (int i=0;i<3;++i) retval.push_back(m_leg1.DV[i]);
	for (int i=0;i<3;++i) retval.push_back(m_leg2.DV[i]);
	return retval;
}

/// Implementation of the sparsity structure.
/**
 * No sparsity present (box-constrained problem).
 */
void sample_return::set_sparsity(int &lenG, std::vector<int> &iGfun, std::vector<int> &jGvar) const
{
	lenG=12;
	iGfun.resize(lenG);
	jGvar.resize(lenG);
	for (int i = 0; i<lenG; ++i)
	{
		iGfun[i] = 0;
		jGvar[i] = i;
	}
}

std::string sample_return::get_name() const
{
	return "Sample return";
}

}} //namespaces

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::sample_return)
