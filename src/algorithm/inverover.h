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

#ifndef PAGMO_ALGORITHM_INVEROVER_H
#define PAGMO_ALGORITHM_INVEROVER_H

#include "../config.h"
#include "../serialization.h"
#include "../population.h"
#include "../problem/tsp.h"
#include "base.h"
#include <algorithm>


namespace pagmo { namespace algorithm {

/// Inver-Over Algorithm (IO)
/**
 * The Inver-Over algorithm is a state-of-the-art genetic algorithm for the Travelling Salesman Problem.
 * It was designed by G. Tao and Z. Michalewicz in 1998.
 * 
 * Note: The algorithm was sightly changed (choice of the next city in a series of inverisons)
 * since with this choice better performance (tour length, computational time) was observed.
 *
 * Note2: The algorithm contains a 2nd stopping criterion depending on the time (generation) an individual was
 * improved. It is possible that more accurate solutions are obtained if this criterion is removed.
 *
 * Note3: The value for the population size is advised to be no smaller than 20.
 * To not have premature convergence, values around 100 are observed to work well.
 *
 * @author Ingmar Getzner (ingmar.getzner@gmail.com)
 */
class __PAGMO_VISIBLE inverover: public base
{
    public:
        inverover(int gen = 500000, double ri = 0.05);

        base_ptr clone() const;
        void evolve(population &) const;
        std::string get_name() const;
	//void set_gen(int gen);
	//int get_gen() const;
	//void set_ri(double ri);
	//double get_ri() const;
    protected:
	size_t compute_idx(const size_t,const size_t,const size_t);

    private:
        friend class boost::serialization::access;
        template <class Archive>
        void serialize(Archive &ar, const unsigned int)
        {
                ar & boost::serialization::base_object<base>(*this);
                ar & m_gen;
                ar & m_ri;
        }

        int m_gen;
        double m_ri;
};

}} //namespaces

BOOST_CLASS_EXPORT_KEY(pagmo::algorithm::inverover)

#endif // PAGMO_ALGORITHM_INVEROVER_H
