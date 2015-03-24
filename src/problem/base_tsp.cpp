/*****************************************************************************
 *   Copyright (C) 2004-2014 The PaGMO development team,                     *
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

#include <algorithm>

#include "base_tsp.h"
#include "../population.h"

namespace pagmo { namespace problem {

    /// Constructor from dimensins and encoding
    /**
     * @param[in] n_cities number of cities
     * @param[in] nc total number of constraints
     * @param[in] nic total number of inequality constraints
     * @param[in] encoding encoding_type, i.e. one of base_tsp::CITIES, base_tsp::FULL, base_tsp::RANDOMKEYS
     */
    base_tsp::base_tsp(int n_cities, int nc, int nic, encoding_type encoding): 
        base(
            (encoding==FULL ? n_cities*(n_cities-1): n_cities), 
            (encoding==RANDOMKEYS ? 0 : (encoding==FULL ? n_cities*(n_cities-1):n_cities)), 
            1, nc, nic, 0.0
        ), 
        m_encoding(encoding), 
        m_n_cities(n_cities)
    {
        switch( m_encoding ) {
            case FULL:
                set_lb(0);
                set_ub(1);
                break;
            case RANDOMKEYS:
                set_lb(0);
                set_ub(1);
                break;
            case CITIES:
                set_lb(0);
                set_ub(m_n_cities-1);
                break;
        }
    }


    /// From FULL to CITIES encoding
    /**
     * Transforms a chromosome in the FULL encoding into a chromosome in the CITIES encoding.
     * If the starting chromosome is unfeasible also the resulting chromosome in the CITIES encoding will be
     * unfeasible.
     *
     * @param[in] x a chromosome in the FULL encoding
     * @return a chromosome in the CITIES encoding
     */
    pagmo::decision_vector base_tsp::full2cities(const pagmo::decision_vector &x) const
    {
        pagmo::decision_vector retval(m_n_cities,0);
        pagmo::population::size_type next_city,cur_city = 0;
        retval[0]=cur_city;
        for(size_t j = 1; j < m_n_cities; j++){
              pagmo::decision_vector::const_iterator iter = std::find(x.begin() + cur_city*(m_n_cities-1), x.begin() + (cur_city+1)*(m_n_cities-1),1);
              next_city = iter - (x.begin() + cur_city*(m_n_cities-1));
              next_city = next_city + ( (next_city >= cur_city) ? 1:0 );
              cur_city=next_city;
              retval.at(j) = std::min(next_city,m_n_cities-1); //the min is to prevent cases where the 1 is not found (unfeasible chromosomes) and thus the city _idx returned would be invalid
        }
        return retval;
    }

    /// From CITIES to FULL encoding
    /**
     * Transforms a chromosome in the CITIES encoding into a chromosome in the FULL encoding.
     * If the starting chromosome is unfeasible also the resulting chromosome in the FULL encoding will be
     * unfeasible.
     *
     * @param[in] x a chromosome in the CITIES encoding
     * @return a chromosome in the FULL encoding
     */
    pagmo::decision_vector base_tsp::cities2full(const pagmo::decision_vector &x) const
    {
        if (x.size() != m_n_cities) 
        {
            pagmo_throw(value_error,"input representation of a tsp solution (CITIES encoding) looks unfeasible [wrong length]");
        }

        pagmo::decision_vector retval(m_n_cities*(m_n_cities-1),0);
        for (std::vector<population::size_type>::size_type i=0; i<x.size()-1; ++i)
        {
            retval.at( x[i]*(m_n_cities-1) + x[i+1] - (x[i+1]>x[i]?1:0) ) = 1;
        } 
        retval.at( x[x.size()-1]*(m_n_cities-1) + x[0] - (x[0]>x[x.size()-1]?1:0) ) = 1;
        return retval;
    }

    bool comparator ( const std::pair<double,int>& l, const std::pair<double,int>& r)
    { return l.first < r.first; }

    /// From RANDOMKEYS to CITIES encoding
    /**
     * Transforms a chromosome in the RANDOMKEYS encoding into a chromosome in the CITIES encoding.
     * If the starting chromosome is unfeasible also the resulting chromosome in the CITIES encoding will be
     * unfeasible.
     *
     * @param[in] x a chromosome in the RANDOMKEYS encoding
     * @return a chromosome in the CITIES encoding
     */

    pagmo::decision_vector base_tsp::randomkeys2cities(const pagmo::decision_vector &x) const
    {
        if (x.size() != m_n_cities) 
        {
            pagmo_throw(value_error,"input representation of a tsp solution (RANDOMKEYS encoding) looks unfeasible [wrong length]");
        }
        pagmo::decision_vector retval(m_n_cities);
        std::vector<std::pair<double,int> > pairs(m_n_cities);
        for (pagmo::decision_vector::size_type i=0;i<m_n_cities;++i) {
            pairs[i].first = x[i];
            pairs[i].second = i;
        }
        std::sort(pairs.begin(),pairs.end(),comparator);
        for (pagmo::decision_vector::size_type i=0;i<m_n_cities;++i) {
            retval[i] = pairs[i].second;
        }
        return retval;
    }

    /// From CITIES to RANDOMKEYS encoding
    /**
     * Transforms a chromosome in the CITIES encoding into a chromosome in the RANDOMKEYS encoding.
     * If the starting chromosome is unfeasible, the resulting chromosome in the RANDOMKEYS encoding will contain zeros,
     * and yet still be feasible. Its inversion using randomkeys2cities will result in a feasible tour.
     *
     * @param[in] x a chromosome in the CITIES encoding
     * @param[in] orig_random_keys a chromosome in the RANDOMKEYS encoding. 
     * @return a chromosome in the RANDOMKEYS encoding
     */
    pagmo::decision_vector base_tsp::cities2randomkeys(const pagmo::decision_vector &cities,const pagmo::decision_vector &orig_random_keys) const
    {
        if (cities.size() != orig_random_keys.size()) 
        {
            pagmo_throw(value_error,"the random keys original vector and the cities vector need to have the same length");
        }
        if (cities.size() != m_n_cities) 
        {
            pagmo_throw(value_error,"input representation of a tsp solution (CITIES encoding) looks unfeasible [wrong length]");
        }
        if ( (*std::max_element(cities.begin(),cities.end()) >= m_n_cities) || (*std::min_element(cities.begin(),cities.end()) < 0) )
        {
            pagmo_throw(value_error,"city indexes outside the allowed bounds");
        }

        pagmo::decision_vector rk(orig_random_keys);
        pagmo::decision_vector retval(m_n_cities);
        std::sort(rk.begin(),rk.end());
        for (pagmo::decision_vector::size_type i=0;i<m_n_cities;++i) {
            retval[cities[i]] = rk[i];
        }
        return retval;
    }

    /// Getter for m_encoding
    /**
     * @return reference to the encoding_type
     */
    base_tsp::encoding_type base_tsp::get_encoding() const
    { 
        return m_encoding; 
    }

    /// Getter for m_n_cities
    /**
     * @return reference to m_n_cities
     */
    decision_vector::size_type base_tsp::get_n_cities() const
    { 
        return m_n_cities; 
    }

}} //namespaces
