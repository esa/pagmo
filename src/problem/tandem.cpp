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
 #include <keplerian_toolbox/epoch.h>

#include "tandem.h"
#include "../AstroToolbox/mga_dsm.h"
#include "../AstroToolbox/misc4Tandem.h"

namespace pagmo { namespace problem {

/// Fly bys considered
const int tandem::Data[24][5] = {
	{3,2,2,2,6},
	{3,2,2,3,6},
	{3,2,2,4,6},
	{3,2,2,5,6},
	{3,2,3,2,6},
	{3,2,3,3,6},
	{3,2,3,4,6},
	{3,2,3,5,6},
	{3,2,4,2,6},
	{3,2,4,3,6},
	{3,2,4,4,6},
	{3,2,4,5,6},
	{3,3,2,2,6},
	{3,3,2,3,6},
	{3,3,2,4,6},
	{3,3,2,5,6},
	{3,3,3,2,6},
	{3,3,3,3,6},
	{3,3,3,4,6},
	{3,3,3,5,6},
	{3,3,4,2,6},
	{3,3,4,3,6},
	{3,3,4,4,6},
	{3,3,4,5,6}
};

/// Dummy sequence (unused)
const int tandem::sequence[5] = {1,1,1,1,1};

/// Constructor
/**
* Instantiates one of the possible TandEM problems
* \param[in] probid This is an integer number from 1 to 24 encoding the fly-by sequence to be used (default is EVEES). Check http://www.esa.int/gsp/ACT/inf/op/globopt/TandEM.htm for more information
* \param[in] tof_ (in years) This is a number setting the constraint on the total time of flight (10 from the GTOP database). If -1 (default) an unconstrained problem is instantiated
*/
tandem::tandem(const int probid, const double tof_):base(18), problem(orbit_insertion,sequence,5,0,0,0,0.98531407996358,80330.0), tof(tof_),copy_of_x(18)
{
	if (probid < 1 || probid > 24) {
		pagmo_throw(value_error,"probid needs to be an integer in [1,24]");
	}
	if (tof_ > 20 && tof_ <5 && tof_!=-1) {
		pagmo_throw(value_error,"time of flight constraint needs to be a number in [5,20] (years)");
	}
	for (int i =0; i<5; i++) {
		problem.sequence[i] = Data[probid-1][i];
	}

	const double lbunc[18] = {5475, 2.50001, 0, 0, 20   , 20  ,  20 , 20  , 0.01, 0.01, 0.01, 0.01, 1.05, 1.05, 1.05, -M_PI, -M_PI, -M_PI};
	const double ubunc[18] = {9132, 4.9, 1, 1, 2500 , 2500, 2500, 2500, 0.99, 0.99, 0.99, 0.99,    10,    10,    10,  M_PI,  M_PI,  M_PI};

	const double lbcon[18] = {5475, 2.50001, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 1.05, 1.05, 1.05, -M_PI, -M_PI, -M_PI};
	const double ubcon[18] = {9132, 4.9, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99,    10,    10,    10,  M_PI,  M_PI,  M_PI};

	if (tof_==-1) set_bounds(lbunc,lbunc+18,ubunc,ubunc+18);
	else set_bounds(lbcon,lbcon+18,ubcon,ubcon+18);
}

/// Clone method.
base_ptr tandem::clone() const
{
	return base_ptr(new tandem(*this));
}

/// Implementation of the objective function.
void tandem::objfun_impl(fitness_vector &f, const decision_vector &x) const
{
	if (tof!=-1) { //constrained problem
		//Here we copy the chromosome into a new vector and we transform its time percentages into days
		copy_of_x = x;
		copy_of_x[4] = x[4]*365.25*tof;
		copy_of_x[5] = x[5]*(365.25*tof-copy_of_x[4]);
		copy_of_x[6] = x[6]*(365.25*tof-copy_of_x[4]-copy_of_x[5]);
		copy_of_x[7] = x[7]*(365.25*tof-copy_of_x[4]-copy_of_x[5]-copy_of_x[6]);
		MGA_DSM(copy_of_x, problem, f[0]);
	} else {	//unconstrained problem
		MGA_DSM(x, problem, f[0]);
	}
	//evaluating the mass from the dvs
	double rE[3];
	double vE[3];
	Planet_Ephemerides_Analytical (x[0],3,rE,vE);
	double VINFE = x[1];
	double udir = x[2];
	double vdir = x[3];
	double vtemp[3];
	vtemp[0]= rE[1]*vE[2]-rE[2]*vE[1];
	vtemp[1]= rE[2]*vE[0]-rE[0]*vE[2];
	vtemp[2]= rE[0]*vE[1]-rE[1]*vE[0];
	double iP1[3];
	double normvE=sqrt(vE[0]*vE[0]+vE[1]*vE[1]+vE[2]*vE[2]);
	iP1[0]=	vE[0]/normvE;
	iP1[1]=	vE[1]/normvE;
	iP1[2]=	vE[2]/normvE;
	double zP1[3];
	double normvtemp=sqrt(vtemp[0]*vtemp[0]+vtemp[1]*vtemp[1]+vtemp[2]*vtemp[2]);
	zP1[0]= vtemp[0]/normvtemp;
	zP1[1]= vtemp[1]/normvtemp;
	zP1[2]= vtemp[2]/normvtemp;
	double jP1[3];
	jP1[0]= zP1[1]*iP1[2]-zP1[2]*iP1[1];
	jP1[1]= zP1[2]*iP1[0]-zP1[0]*iP1[2];
	jP1[2]= zP1[0]*iP1[1]-zP1[1]*iP1[0];
	double theta=2*M_PI*udir; 		//See Picking a Point on a Sphere
	double phi=acos(2*vdir-1)-M_PI/2; //In this way: -pi/2<phi<pi/2 so phi can be used as out-of-plane rotation
	double vinf[3];
	vinf[0]=VINFE*(cos(theta)*cos(phi)*iP1[0]+sin(theta)*cos(phi)*jP1[0]+sin(phi)*zP1[0]);
	vinf[1]=VINFE*(cos(theta)*cos(phi)*iP1[1]+sin(theta)*cos(phi)*jP1[1]+sin(phi)*zP1[1]);
	vinf[2]=VINFE*(cos(theta)*cos(phi)*iP1[2]+sin(theta)*cos(phi)*jP1[2]+sin(phi)*zP1[2]);
	//We rotate it to the equatorial plane
	ecl2equ(vinf,vinf);
	//And we find the declination in degrees
	double normvinf=sqrt(vinf[0]*vinf[0]+vinf[1]*vinf[1]+vinf[2]*vinf[2]);
	double sindelta = vinf[2] / normvinf;
	double declination = asin(sindelta)/M_PI*180;

	//double m_initial = SoyuzFregat(VINFE,declination);
	double m_initial = Atlas501(VINFE,declination);

	//We evaluate the final mass
	double Isp = 312;
	double g0 = 9.80665;
	double sumDVvec=0;

	for(unsigned int i=1;i<=5;i++) {
		sumDVvec=sumDVvec+problem.DV[i];
	}
	double m_final;
	sumDVvec=sumDVvec+0.165; //losses for 3 swgbys + insertion
	m_final = m_initial * exp(-sumDVvec/Isp/g0*1000);
	f[0] = -log(m_final);
}

/// Outputs a stream of the trajectory data
/**
 * While the chromosome contains all necessary information to describe a trajectory, mission analysits
 * often require a different set of data to understand a trajectory. This method outputs a stream with
 * information on the trajectory that is otherwise 'hidden' in the chromosome
 *
 * \param[in] x chromosome representing the trajectory in the optimization process
 * \returns an std::string with launch dates, DV magnitues and other information on the trajectory
 */
std::string tandem::pretty(const std::vector<double> &x) const
{
	double obj=0;
	std::vector<double> printablex;
	if (tof!=-1) { //constrained problem
		//Here we copy the chromosome into a new vector and we transform its time percentages into days
		copy_of_x = x;
		copy_of_x[4] = x[4]*365.25*tof;
		copy_of_x[5] = x[5]*(365.25*tof-copy_of_x[4]);
		copy_of_x[6] = x[6]*(365.25*tof-copy_of_x[4]-copy_of_x[5]);
		copy_of_x[7] = x[7]*(365.25*tof-copy_of_x[4]-copy_of_x[5]-copy_of_x[6]);
		MGA_DSM(copy_of_x, problem, obj);
		printablex=copy_of_x;
	} else {	//unconstrained problem
		MGA_DSM(x, problem, obj);
		printablex=x;
	}
	std::ostringstream s;
	s.precision(15);
	s << std::scientific;
	const size_t seq_size = (x.size() + 2) / 4;
	pagmo_assert((x.size() + 2) % 4 == 0 && seq_size >= 2);
	pagmo_assert(problem.sequence.size() == seq_size);
	s << "Flyby sequence:\t\t\t";
	for (size_t i = 0; i < seq_size; ++i) {
		s << problem.sequence[i];
	}
	s << '\n';
	s << "Departure epoch (mjd2000):\t" << printablex[0] << '\n';
	s << "Departure epoch:\t\t" << ::kep_toolbox::epoch(printablex[0],::kep_toolbox::epoch::MJD2000) << '\n';
	s << "Vinf polar components:\t\t";
	for (size_t i = 0; i < 3; ++i) {
		s << printablex[i + 1] << ' ';
	}
	s << '\n';
	double totaltime = 0;
	for (size_t i = 0; i < seq_size - 1; ++i) {
		s << "Leg time of flight:\t\t" << printablex[i + 4] << '\n';
		totaltime += printablex[i + 4];
	}
	s << "Total time of flight:\t\t" << totaltime << '\n';
	for (size_t i = 0; i < seq_size - 2; ++i) {
		s << "Flyby radius:\t\t\t" << printablex[i + 2 * (seq_size + 1)] << '\n';
	}
	totaltime=printablex[0];
	for (size_t i = 0; i < seq_size - 2; ++i) {
	totaltime += printablex[i + 4];
		s << "Flyby date:\t\t\t" << ::kep_toolbox::epoch(totaltime,::kep_toolbox::epoch::MJD2000) << '\n';
	}
	for (size_t i = 0; i < seq_size - 2; ++i) {
		s << "Vinf at flyby:\t\t\t" << std::sqrt(problem.vrelin_vec[i]) << '\n';
	}
	for (size_t i = 0; i < seq_size - 1; ++i) {
		s << "dsm" << i+1 << ":\t\t\t\t" << problem.DV[i+1] << '\n';
	}
	s << "Final DV:\t\t\t" << problem.DV.back() << '\n';
	return s.str();
}

/// Implementation of the sparsity structure.
/**
 * This implementation does not make use of problem::base::estimate_sparsity
 * as the numerical difficulties introduced by the objective function definition through a logarithm
 * makes automated detection unreliable (e.g. also SNOPT algorithm fails).
 * Clearly, as the problem is box constrained no sarsity is present.
 */

void tandem::set_sparsity(int &lenG, std::vector<int> &iGfun, std::vector<int> &jGvar) const
{
	lenG=18;
	iGfun.resize(18);
	jGvar.resize(18);
	for (int i = 0; i<lenG; ++i)
	{
		iGfun[i] = 0;
		jGvar[i] = i;
	}
}

std::string tandem::get_name() const
{
	return "TandEM";
}

}} //namespaces

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::tandem)
