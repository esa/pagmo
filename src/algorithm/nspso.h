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

#ifndef PAGMO_ALGORITHM_NSPSO_H
#define PAGMO_ALGORITHM_NSPSO_H

#include "../config.h"
#include "../serialization.h"
#include "base.h"

namespace pagmo { namespace algorithm {

/// Non-dominated Sorting Particle Swarm Optimizer (NSPSO)
/**
 *
 * asd
 *
 * @author Andrea Mambrini (andrea.mambrini@gmail.com)
 * @author Annalisa Riccardi (nina1983@gmail.com)
 *
 * @see "Xiaodong Li - A Non-dominated Sorting Particle Swarm Optimizer for Multiobjective Optimization"
 **/

class __PAGMO_VISIBLE nspso: public base
{
public:
	nspso(int gen=10, double minW = 0.4, double maxW = 1.0, double C1 = 2.0, double C2 = 2.0,
		  double CHI = 1.0, double v_coeff = 0.5, int leader_selection_range = 30);
	nspso(const nspso &);

	base_ptr clone() const;
	void evolve(population &) const;
	std::string get_name() const;

protected:
	std::string human_readable_extra() const;

private:
	void compute_niche_count(std::vector<int> &, const std::vector<std::vector<double> > &, double) const;
	double euclidian_distance(const std::vector<double> &, const std::vector<double> &) const;
	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive &ar, const unsigned int)
	{
		ar & boost::serialization::base_object<base>(*this);
		ar & const_cast<int &>(m_gen);
		ar & const_cast<double &>(m_minW);
		ar & const_cast<double &>(m_maxW);
		ar & const_cast<double &>(m_C1);
		ar & const_cast<double &>(m_C2);
		ar & const_cast<double &>(m_CHI);
		ar & const_cast<double &>(m_v_coeff);
		ar & const_cast<int &>(m_leader_selection_range);

	}
	//Number of generations
	const int m_gen;
	const double m_minW;
	const double m_maxW;
	const double m_C1;
	const double m_C2;
	const double m_CHI;
	const double m_v_coeff;
	const int m_leader_selection_range;

};

}} //namespaces

BOOST_CLASS_EXPORT_KEY(pagmo::algorithm::nspso);

#endif // PAGMO_ALGORITHM_NSPSO_H
