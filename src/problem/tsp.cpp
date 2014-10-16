/*****************************************************************************
 *   Copyright (C) 2004-2014 The PaGMO development team,                     *
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

#include "tsp.h"
#include "../population.h"

namespace pagmo { namespace problem {

    /// Default constructor
    /**
     * This constructs a 3-cities symmetric problem (naive TSP) 
     * with weight matrix [[0,1,1][1,0,1][1,1,0]] and FULL encoding
     */
    tsp::tsp(): base(6, 6, 1, 8, 2, 0.0), m_n_cities(3), m_weights(), m_encoding(tsp::FULL)
    {
        std::vector<double> dumb(3,0);
        m_weights = std::vector<std::vector<double> > (3,dumb);
        m_weights[0][1] = 1;
        m_weights[0][2] = 1;
        m_weights[2][1] = 1;
        m_weights[1][0] = 1;
        m_weights[2][0] = 1;
        m_weights[1][2] = 1;
        set_lb(0);
        set_ub(1);
    }

    /// Constructor from weight matrix and encoding
    /**
     * Constructs a TSP with the input weight matrix and the selected encoding
     * @param[in] weights an std::vector of std::vector representing a square matrix.
     * @param[in] encoding_ a pagmo::problem::tsp::encoding representing the chosen encoding
     */
    tsp::tsp(const std::vector<std::vector<double> >& weights, const encoding& encoding_): 
        base(
            tsp::compute_dimensions(weights[0].size(),encoding_)[0],
            tsp::compute_dimensions(weights[0].size(),encoding_)[1],
            tsp::compute_dimensions(weights[0].size(),encoding_)[2],
            tsp::compute_dimensions(weights[0].size(),encoding_)[3],
            tsp::compute_dimensions(weights[0].size(),encoding_)[4],
            (double)tsp::compute_dimensions(weights[0].size(),encoding_)[5]
        ), m_n_cities(weights[0].size()), m_weights(weights), m_encoding(encoding_)
    {
        check_weights(m_weights);
        switch( encoding_ ) {
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

    /// Clone method.
    base_ptr tsp::clone() const
    {
        return base_ptr(new tsp(*this));
    }

    /// Checks if we can instantiate a TSP or ATSP problem
    /**
     * Checks if a matrix (std::vector<std::vector<double>>) 
     * is square or bidirectional (e.g. no one way links between vertices).
     * If none of the two conditions are true, we can not have a tsp problem.
     * @param matrix - the adjacency matrix (two dimensional std::vector)
     * @throws pagmo_throw - matrix is not square and/or graph is not bidirectional
     */
    void tsp::check_weights(const std::vector<std::vector<double> > &matrix) const 
    {   
        decision_vector::size_type n_cols = matrix.size();
        
        for (decision_vector::size_type i = 0; i < n_cols; ++i) {
            decision_vector::size_type n_rows = matrix.at(i).size();
            // check if the matrix is square
            if (n_rows != n_cols)
                pagmo_throw(value_error, "adjacency matrix is not square");
            
            for (size_t j = 0; j < n_rows; ++j) {
                if (i == j && matrix.at(i).at(j) != 0)
                    pagmo_throw(value_error, "main diagonal elements must all be zeros.");
                if (i != j && !matrix.at(i).at(j)) // fully connected
                    pagmo_throw(value_error, "adjacency matrix contains zero values.");
                if (i != j && (!matrix.at(i).at(j)) == matrix.at(i).at(j)) // fully connected
                    pagmo_throw(value_error, "adjacency matrix contains NaN values.");                    
            }
        }
    }

    boost::array<int, 6> tsp::compute_dimensions(decision_vector::size_type n_cities, encoding encoding_)
    {
        boost::array<int,6> retval;
        switch( encoding_ ) {
            case FULL:
                retval[0] = n_cities*(n_cities-1);
                retval[1] = n_cities*(n_cities-1);
                retval[2] = 1;
                retval[3] = n_cities*(n_cities-1)+2;
                retval[4] = (n_cities-1)*(n_cities-2);
                retval[5] = 0.0;
                break;
            case RANDOMKEYS:
                retval[0] = n_cities;
                retval[1] = 0;
                retval[2] = 1;
                retval[3] = 0;
                retval[4] = 0;
                retval[5] = 0.0;
                break;
            case CITIES:
                retval[0] = n_cities;
                retval[1] = n_cities;
                retval[2] = 1;
                retval[3] = 1;
                retval[4] = 0;
                retval[5] = 0.0;
                break;
        }
        return retval;
    }

    void tsp::objfun_impl(fitness_vector &f, const decision_vector& x) const 
    {
        decision_vector tour;
        f[0]=0;
        switch( m_encoding ) {
            case FULL:
                tour = full2cities(x);
                break;
            case RANDOMKEYS:
                tour = randomkeys2cities(x);
                break;
            case CITIES:
                tour = x;
                break;
        }
        for (decision_vector::size_type i=0; i<m_n_cities-1; ++i) {
            f[0] += m_weights.at(tour.at(i)).at(tour.at(i+1));
        }
        f[0]+= m_weights.at(tour.at(m_n_cities-1)).at(tour.at(0));
    }

    size_t compute_idx(const size_t i, const size_t j, const size_t n) 
    {
        pagmo_assert( i!=j && i<n && j<n );
        return i*(n-1) + j - (j>i? 1:0);
    }

    void tsp::compute_constraints_impl(constraint_vector &c, const decision_vector& x) const 
    {
        decision_vector::size_type n = get_n_cities();

        switch( m_encoding ) 
        {
            case FULL:
            {
                // 1 - We set the equality constraints
                for (size_t i = 0; i < n; i++) {
                    c[i] = 0;
                    c[i+n] = 0;
                    for (size_t j = 0; j < n; j++) {
                        if(i==j) continue; // ignoring main diagonal
                        decision_vector::size_type rows = compute_idx(i, j, n);
                        decision_vector::size_type cols = compute_idx(j, i, n);
                        c[i] += x[rows];
                        c[i+n] += x[cols];
                    }
                    c[i] = c[i]-1;
                    c[i+n] = c[i+n]-1;
                }

                //2 - We set the inequality constraints
                //2.1 - First we compute the uj (see http://en.wikipedia.org/wiki/Travelling_salesman_problem#Integer_linear_programming_formulation)
                //      we start always out tour from the first city, without loosing generality
                size_t next_city = 0,current_city = 0;
                std::vector<int> u(n);
                for (size_t i = 0; i < n; i++) {
                    u[current_city] = i+1;
                    for (size_t j = 0; j < n; j++) 
                    {
                        if (current_city==j) continue;
                        if (x[compute_idx(current_city, j, n)] == 1) 
                        {
                            next_city = j;
                            break;
                        }
                    }
                    current_city = next_city;
                }
                int count=0;
                for (size_t i = 1; i < n; i++) {
                    for (size_t j = 1; j < n; j++) 
                    {
                        if (i==j) continue;
                        c[2*n+count] = u[i]-u[j] + (n+1) * x[compute_idx(i, j, n)] - n;
                        count++;
                    }
                }
                break;
            }
            case RANDOMKEYS:
                break;
            case CITIES:
            {
                std::vector<population::size_type> range(m_n_cities);
                for (std::vector<population::size_type>::size_type i=0; i<range.size(); ++i) 
                {
                    range[i]=i;
                }
                c[0] = !std::is_permutation(x.begin(),x.end(),range.begin());
                break;
            }
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
    pagmo::decision_vector tsp::full2cities(const pagmo::decision_vector &x) const
    {
        if (x.size() != (m_n_cities-1)*m_n_cities )
        {
            pagmo_throw(value_error,"input representation of a tsp solution (FULL encoding) has the wrong length");
        }

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
    pagmo::decision_vector tsp::cities2full(const pagmo::decision_vector &x) const
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

    pagmo::decision_vector tsp::randomkeys2cities(const pagmo::decision_vector &x) const
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
    pagmo::decision_vector tsp::cities2randomkeys(const pagmo::decision_vector &cities,const pagmo::decision_vector &orig_random_keys) const
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

    /// Getter for m_weights
    /**
     * @return reference to m_weights
     */
    const std::vector<std::vector<double> >&  tsp::get_weights() const
    { 
        return m_weights; 
    }

    /// Getter for m_encoding
    /**
     * @return reference to m_encoding
     */
    const tsp::encoding& tsp::get_encoding() const
    { 
        return m_encoding; 
    }
    
    /// Getter for m_n_cities
    /**
     * @return reference to m_n_cities
     */
    const decision_vector::size_type& tsp::get_n_cities() const
    { 
        return m_n_cities; 
    }

    /// Returns the problem name
    std::string tsp::get_name() const
    {
        return "Travelling Salesman Problem (TSP-ATSP)";
    }

    /// Extra human readable info for the problem.
    /**
     * @return a std::string containing a list of vertices and edges
     */
    std::string tsp::human_readable_extra() const 
    {
        std::ostringstream oss;
        oss << "\n\tNumber of cities: " << m_n_cities << '\n';
        oss << "\tEncoding: ";
        switch( m_encoding ) {
            case FULL:
                oss << "FULL" << '\n';
                break;
            case RANDOMKEYS:
                oss << "RANDOMKEYS" << '\n';
                break;
            case CITIES:
                oss << "CITIES" << '\n';
                break;
        }
        oss << "\tWeight Matrix: \n";
        for (decision_vector::size_type i=0; i<m_n_cities; ++i)
        {
            oss << "\t\t" << m_weights[i] << '\n';
            if (i>5)
            {
                oss << "\t\t..." << '\n';
                break;
            }
        }
        return oss.str();
    }

    
}} //namespaces