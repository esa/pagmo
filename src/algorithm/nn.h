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

#ifndef PAGMO_ALGORITHM_IO_H
#define PAGMO_ALGORITHM_IO_H

#include "../config.h"
#include "../serialization.h"
#include "../population.h"
#include "../problem/tsp.h"
#include "base.h"
#include <algorithm>


namespace pagmo { namespace algorithm {

/// Nearest Neighbor Algorithm (NN)
/**
 * The Nearest Neighbor algorithm generates a tour starting either from a single, in the input, specified vertex
 * or loops over all possible initial vertices, computes the corresponding tours and returns the shortest tour.
 *
 * @author Ingmar Getzner (ingmar.getzner@gmail.com)
 */
class __PAGMO_VISIBLE nn: public base
{
    public:
        nn(int start_city = -1);

        base_ptr clone() const;
        void evolve(population &) const;
        std::string get_name() const;

    private:
        friend class boost::serialization::access;
        template <class Archive>
        void serialize(Archive &ar, const unsigned int)
        {
                ar & boost::serialization::base_object<base>(*this);
                ar & m_start_city;
        }

        int m_start_city;
};

}} //namespaces

BOOST_CLASS_EXPORT_KEY(pagmo::algorithm::nn)

#endif // PAGMO_ALGORITHM_NN_H
