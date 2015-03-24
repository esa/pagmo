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

#include "tsp_vrplc.h"
#include "../population.h"

namespace pagmo { namespace problem {

    /// Default constructor
    /**
     * This constructs a 3-cities symmetric problem (naive TSP) 
     * with weight matrix [[0,1,1][1,0,1][1,1,0]] and RANDOMKEYS encoding
     */
    tsp_vrplc::tsp_vrplc() : base_tsp(3, 0, 0 , base_tsp::RANDOMKEYS), m_weights(), m_capacity(1.1)
    {
        std::vector<double> dumb(3,0);
        m_weights = std::vector<std::vector<double> > (3,dumb);
        m_weights[0][1] = 1;
        m_weights[0][2] = 1;
        m_weights[2][1] = 1;
        m_weights[1][0] = 1;
        m_weights[2][0] = 1;
        m_weights[1][2] = 1;
    }

    /// Constructor from weight matrix, encoding and capacity
    /**
     * Constructs a TSP with the input weight matrix and the selected encoding
     * @param[in] weights an std::vector of std::vector representing a square matrix.
     * @param[in] encoding a pagmo::problem::tsp::encoding representing the chosen encoding
     * @param[in] capacity maximum vehicle capacity
     */
    tsp_vrplc::tsp_vrplc(const std::vector<std::vector<double> >& weights, const base_tsp::encoding_type& encoding, const double& capacity): 
        base_tsp(weights.size(), 
            compute_dimensions(weights.size(), encoding)[0],
            compute_dimensions(weights.size(), encoding)[1],
            encoding
        ),  m_weights(weights), m_capacity(capacity)
    {
        if (m_capacity <= 0)
        {
            pagmo_throw(value_error, "Maximum vehicle capacity needs to be strictly positive");
        }
        check_weights(m_weights);
    }

    /// Clone method.
    base_ptr tsp_vrplc::clone() const
    {
        return base_ptr(new tsp_vrplc(*this));
    }

    /// Checks if we can instantiate a TSP or ATSP problem
    /**
     * Checks if a matrix (std::vector<std::vector<double>>) 
     * is square or bidirectional (e.g. no one way links between vertices).
     * If none of the two conditions are true, we can not have a tsp problem.
     *
     * @param matrix the adjacency matrix (two dimensional std::vector)
     * @throws pagmo_throw matrix is not square and/or graph is not bidirectional
     */
    void tsp_vrplc::check_weights(const std::vector<std::vector<double> > &matrix) const 
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

    boost::array<int, 2> tsp_vrplc::compute_dimensions(decision_vector::size_type n_cities, base_tsp::encoding_type encoding)
    {
        boost::array<int,2> retval;
        switch( encoding ) {
            case FULL:
                retval[0] = n_cities*(n_cities-1)+2;
                retval[1] = (n_cities-1)*(n_cities-2);
                break;
            case RANDOMKEYS:
                retval[0] = 0;
                retval[1] = 0;
                break;
            case CITIES:
                retval[0] = 1;
                retval[1] = 0;
                break;
        }
        return retval;
    }

    void tsp_vrplc::objfun_impl(fitness_vector &f, const decision_vector& x) const 
    {
        f[0] = 1;
        double stl = 0;
        decision_vector tour;
        decision_vector::size_type n_cities = get_n_cities();
        switch( get_encoding() ) {
            case FULL:
            {
                tour = full2cities(x);
                break;
            }
            case RANDOMKEYS:
            {
                tour = randomkeys2cities(x);
                break;
            }
            case CITIES:
            {
                tour = x;
                break;
            }
        }
        for (decision_vector::size_type i=0; i<n_cities-1; ++i) {
            stl += m_weights[tour[i]][tour[i+1]];
            if(stl > m_capacity)
            {
                stl = 0;
                f[0] += 1;
            }
            else
            {
                f[0] += (m_weights[tour[i]][tour[i+1]])/(n_cities*m_capacity);
            }
        }
        return;
    }
    
    /// Returns the tours
    /**
     * This function takes an Hamiltonian path and brakes it down in the sub-tours
     *
     * @param[in] x a tour in the CITIES encoding
     * @return an std::vector<std::vector<double> > containing the subtours (city indexes)
     */
    std::vector<std::vector<double> > tsp_vrplc::return_tours(const decision_vector& x) const
    {
        std::vector<std::vector<double> > tours;
        double stl = 0;
        decision_vector::size_type n_cities = get_n_cities();
        std::vector<double> cur_tour;
        for (decision_vector::size_type i=0; i<n_cities-1; ++i) 
        {
            cur_tour.push_back(x[i]);
            stl += m_weights[x[i]][x[i+1]];
            if(stl > m_capacity)
            {
                    stl = 0;
                    tours.push_back(cur_tour);
                    cur_tour.erase(cur_tour.begin(),cur_tour.end());
            }
        }
        cur_tour.push_back(x[n_cities-1]);
        tours.push_back(cur_tour);
        return tours;
    }
    

    size_t tsp_vrplc::compute_idx(const size_t i, const size_t j, const size_t n) const
    {
        pagmo_assert( i!=j && i<n && j<n );
        return i*(n-1) + j - (j>i? 1:0);
    }

    void tsp_vrplc::compute_constraints_impl(constraint_vector &c, const decision_vector& x) const 
    {
        decision_vector::size_type n_cities = get_n_cities();

        switch( get_encoding() ) 
        {
            case FULL:
            {
                // 1 - We set the equality constraints
                for (size_t i = 0; i < n_cities; i++) {
                    c[i] = 0;
                    c[i+n_cities] = 0;
                    for (size_t j = 0; j < n_cities; j++) 
                    {
                        if(i==j) continue; // ignoring main diagonal
                        decision_vector::size_type rows = compute_idx(i, j, n_cities);
                        decision_vector::size_type cols = compute_idx(j, i, n_cities);
                        c[i] += x[rows];
                        c[i+n_cities] += x[cols];
                    }
                    c[i] = c[i]-1;
                    c[i+n_cities] = c[i+n_cities]-1;
                }

                //2 - We set the inequality constraints
                //2.1 - First we compute the uj (see http://en.wikipedia.org/wiki/Travelling_salesman_problem#Integer_linear_programming_formulation)
                //      we start always out tour from the first city, without loosing generality
                size_t next_city = 0,current_city = 0;
                std::vector<int> u(n_cities);
                for (size_t i = 0; i < n_cities; i++) 
                {
                    u[current_city] = i+1;
                    for (size_t j = 0; j < n_cities; j++) 
                    {
                        if (current_city==j) continue;
                        if (x[compute_idx(current_city, j, n_cities)] == 1) 
                        {
                            next_city = j;
                            break;
                        }
                    }
                    current_city = next_city;
                }
                int count=0;
                for (size_t i = 1; i < n_cities; i++) {
                    for (size_t j = 1; j < n_cities; j++) 
                    {
                        if (i==j) continue;
                        c[2*n_cities+count] = u[i]-u[j] + (n_cities+1) * x[compute_idx(i, j, n_cities)] - n_cities;
                        count++;
                    }
                }
                break;
            }
            case RANDOMKEYS:
                break;
            case CITIES:
            {
                std::vector<population::size_type> range(n_cities);
                for (std::vector<population::size_type>::size_type i=0; i<range.size(); ++i) 
                {
                    range[i]=i;
                }
                c[0] = !std::is_permutation(x.begin(),x.end(),range.begin());
                break;
            }
        }
        return;
    }

    /// Definition of the distance function
    double tsp_vrplc::distance(decision_vector::size_type i, decision_vector::size_type j) const
    {
        return m_weights[i][j];
    }

    /// Getter for m_weights
    /**
     * @return reference to m_weights
     */
    const std::vector<std::vector<double> >&  tsp_vrplc::get_weights() const
    { 
        return m_weights; 
    }

    /// Getter for m_capacity
    /**
     * @return reference to m_capacity
     */
    const double&  tsp_vrplc::get_capacity() const
    { 
        return m_capacity; 
    }

    /// Returns the problem name
    std::string tsp_vrplc::get_name() const
    {
        return "Vehicle Routing Problem with Limited Vehicle Capacity";
    }

    /// Extra human readable info for the problem.
    /**
     * @return a std::string containing a list of vertices and edges
     */
    std::string tsp_vrplc::human_readable_extra() const 
    {
        std::ostringstream oss;
        oss << "\n\tNumber of cities: " << get_n_cities() << '\n';
        oss << "\tEncoding: ";
        switch( get_encoding()  ) {
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
        oss << "\tMaximum vehicle capacity: " << m_capacity << std::endl;
        oss << "\tWeight Matrix: \n";
        for (decision_vector::size_type i=0; i<get_n_cities() ; ++i)
        {
            oss << "\t\t" << m_weights.at(i) << '\n';
            if (i>5)
            {
                oss << "\t\t..." << '\n';
                break;
            }
        }
        return oss.str();
    }

    
}} //namespaces

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::tsp_vrplc)
