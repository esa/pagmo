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

#ifndef PAGMO_ALGORITHM_SMS_EMOA_H
#define PAGMO_ALGORITHM_SMS_EMOA_H

#include "base.h"
#include "../config.h"
#include "../serialization.h"
#include "../util/hypervolume.h"
#include "../util/hv_algorithm/base.h"


namespace pagmo { namespace algorithm {

/// S-metric selection evolutionary multiobjective optimisation algorithm (SMS-EMOA)
/**
 * SMS-EMOA is a S-metric (hypervolume indicator) based evolutionary algorithm.
 *
 * @see Nicola Beume, Boris Naujoks, Michael Emmerich, "SMS-EMOA: Multiobjective selection based on dominated hypervolume"
 *
 * @author Krzysztof Nowak kn@kiryx.net
 */

class __PAGMO_VISIBLE sms_emoa: public base
{
public:
	sms_emoa(const sms_emoa &);
	sms_emoa(int gen=100, int sel_m = 2, double cr = 0.95, double eta_c = 10, double m = 0.01, double eta_m = 50);
	sms_emoa(pagmo::util::hv_algorithm::base_ptr hv_algorithm, int gen=100, int sel_m = 2, double cr = 0.95, double eta_c = 10, double m = 0.01, double eta_m = 50);
	base_ptr clone() const;
	pagmo::util::hv_algorithm::base_ptr get_hv_algorithm() const;
	void evolve(population &) const;
	std::string get_name() const;
	
protected:
	std::string human_readable_extra() const;
	
private:
	void validate_parameters();
	void crossover(decision_vector&, decision_vector&, pagmo::population::size_type, pagmo::population::size_type,const pagmo::population&) const;
	void mutate(decision_vector&, const pagmo::population&) const;
	population::size_type evaluate_s_metric_selection(const population & pop) const;
	
	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive &ar, const unsigned int)
	{
		ar & boost::serialization::base_object<base>(*this);
		ar & m_hv_algorithm;
		ar & const_cast<int &>(m_gen);
		ar & const_cast<int &>(m_sel_m);
		ar & const_cast<double &>(m_cr);
		ar & const_cast<double &>(m_eta_c);
		ar & const_cast<double &>(m_m);
		ar & const_cast<double &>(m_eta_m);
	}
	// Algorithm used for computation of the least contributor
	pagmo::util::hv_algorithm::base_ptr m_hv_algorithm;
	//Number of generations
	const int m_gen;
	// Selection method for dominated fronts
	// 1 - always use least contributor for every front
	// 2 - use domination count for points in fronts other than front[0]
	const int m_sel_m;
	//Crossover rate
	const double m_cr;
	// Ditribution index for crossover
	const double m_eta_c;
	// Mutation rate
	const double m_m;
	// Ditribution index for mutation
	const double m_eta_m;

};


}} //namespaces

BOOST_CLASS_EXPORT_KEY(pagmo::algorithm::sms_emoa)

#endif // PAGMO_ALGORITHM_SMS_EMOA_H
