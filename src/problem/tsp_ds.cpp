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

#include <numeric>
#include <keplerian_toolbox/astro_constants.h>
#include <keplerian_toolbox/lambert_problem.h>
#include <keplerian_toolbox/core_functions/array3D_operations.h>
#include <keplerian_toolbox/core_functions/ic2par.h>
#include <keplerian_toolbox/epoch.h>

#include "tsp_ds.h"
#include "../exceptions.h"
#include "../population.h"



namespace pagmo { namespace problem {

    /// Constructor
    /**
     * Constructs a TSP - Debris Selection Problem 
     *
     * @param[in] planets         an std::vector of kep_toolbox::planet_ptr containing the orbiting objects one wants to visit
     * @param[in] values          a std::vector representing the value (score) of a visit (is the same order as planets)
     * @param[in] max_DV          the maximum DV available on board the spacecraft (corresponding to max_path_length parameter in the problem::tsp_cs problem)
     * @param[in] times           an std::vector of doubles containing [t0,T1,T2,...] the starting mission epoch and the time of flights between orbiting objects
     * @param[in] encoding        a pagmo::problem::tsp::encoding representing the chosen encoding
     */
    tsp_ds::tsp_ds (            
        const std::vector<kep_toolbox::planet::planet_ptr>& planets, 
        const std::vector<double>& values,
        const double max_DV, 
        const std::vector<double>&  epochs, 
        const base_tsp::encoding_type & encoding
        ) : base_tsp(
                planets.size(), 
                compute_dimensions(planets.size(), encoding)[0],
                compute_dimensions(planets.size(), encoding)[1],
                encoding
        ),  m_planets(planets), m_values(values), m_max_DV(max_DV), m_epochs(epochs), m_mu(planets[0]->get_mu_central_body()), m_DV(values.size() - 1)
    {
        if (planets.size() != values.size()){
            pagmo_throw(value_error, "Planet list must have the same size as values list");
        }
        if (planets.size() != epochs.size()){
            pagmo_throw(value_error, "Planet list must have the same size as times list");
        }
        if (planets.size() < 3 ){
            pagmo_throw(value_error, "Planet list must contain at least 3 elements");
        }

        for (auto pl_ptr : m_planets) 
        {
            if (pl_ptr->get_mu_central_body() != m_mu) 
            {
               pagmo_throw(value_error, "All planets in planet list must have the same gravity constant"); 
            }
        }
        precompute_ephemerides();
    }

    /// Clone method.
    base_ptr tsp_ds::clone() const
    {
        return base_ptr(new tsp_ds(*this));
    }

    void tsp_ds::objfun_impl(fitness_vector &f, const decision_vector& x) const 
    {
        f[0]=0;
        decision_vector tour;
        decision_vector::size_type n_cities = get_n_cities();
        double cum_p, ham_path_len;
        decision_vector::size_type dumb1, dumb2;

        switch( get_encoding() ) {
            case FULL:
            {
                tour = full2cities(x);
                find_subsequence(tour, cum_p, ham_path_len, dumb1, dumb2);
                break;
            }
            case RANDOMKEYS:
            {
                tour = randomkeys2cities(x);
                find_subsequence(tour, cum_p, ham_path_len, dumb1, dumb2);
                break;
           }
            case CITIES:
           {
                find_subsequence(x, cum_p, ham_path_len, dumb1, dumb2);
           }
        }
        ham_path_len = std::accumulate(m_DV.begin(), m_DV.end(), 0.0);
        f[0] = -(cum_p + 1 - ham_path_len / (n_cities * m_max_DV * 10));    
        return;
    }

    /// Finds, in a hamiltonian path, the best subsequence satisfying the max_DV constraint
    /**
     * Finds, in a hamiltonian path, the best subsequence satisfying the max_DV constraint
     * If the input planet sequence does not represent an Hamiltonian path, (i.e. its an unfeasible chromosome)
     * the algorithm behaviour is undefined
     *
     * @param[in]  tour the hamiltonian path encoded with a CITIES encoding (i.e. the list of cities ids)
     * @param[out] retval_p the total cumulative value of the subpath
     * @param[out] retval_l the total saved length of the Hamiltonian path
     * @param[out] retval_it_l the id of the city where the subpath starts
     * @param[out] retval_it_r the id of the city where the subpath ends
     * @throws value_error if the input tour length is not equal to the city number

     */
    void tsp_ds::find_subsequence(const decision_vector& tour, double& retval_p, double& retval_l, decision_vector::size_type& retval_it_l, decision_vector::size_type& retval_it_r, const bool static_computations) const
    {
        if (tour.size() != get_n_cities())
        {
            pagmo_throw(value_error, "tour dimension must be equal to the number of planets");
        }

        // We declare the necessary variable
        decision_vector::size_type n_cities = get_n_cities();
        decision_vector::size_type it_l = 0, it_r = 0;
        bool cond_r = true, cond_l = true;
        double cum_p = m_values[tour[0]];
        double saved_length = m_max_DV;

        // We initialize the starting values
        retval_p = cum_p;
        retval_l = saved_length;
        retval_it_l = it_l;
        retval_it_r = it_r;

        // We precompute all DVs (stored in m_DV)
        compute_DVs(tour, static_computations);

        // Main body of the double loop
        while(cond_l)
        {
            while(cond_r) 
            {
                // We increment the right "pointer" updating the value and length of the path
                saved_length -= m_DV[it_r];
                cum_p += m_values[tour[(it_r + 1)]];
                it_r += 1;

                // We update the various retvals only if the new subpath is valid
                if (saved_length < 0 || (it_l == it_r))
                {
                    cond_r = false;
                }
                else if (cum_p > retval_p)
                {
                    retval_p = cum_p;
                    retval_l = saved_length;
                    retval_it_l = it_l;
                    retval_it_r = it_r;
                }
                else if (cum_p == retval_p)
                {
                    if (saved_length > retval_l)
                    {
                        retval_p = cum_p;
                        retval_l = saved_length;
                        retval_it_l = it_l;
                        retval_it_r = it_r;
                    }
                }

                if (it_r==n_cities-1) 
                {
                    goto EndOfLoop;
                }
            }
            // We get out if all cities are included in the current path
            if (it_l == it_r)
            {
                cond_l = false;
            }
            else
            {
                // We increment the left "pointer" updating the value and length of the path
                saved_length += m_DV[it_l];
                cum_p -= m_values[tour[it_l]];
                it_l += 1;
                // We update the various retvals only if the new subpath is valid
                if (saved_length > 0)
                {
                    cond_r = true;
                    if (cum_p > retval_p)
                    {
                        retval_p = cum_p;
                        retval_l = saved_length;
                        retval_it_l = it_l;
                        retval_it_r = it_r;
                    }
                    else if (cum_p == retval_p)
                    {
                        if (saved_length > retval_l)
                        {
                            retval_p = cum_p;
                            retval_l = saved_length;
                            retval_it_l = it_l;
                            retval_it_r = it_r;
                        }
                    }
                }
                if (it_l == n_cities)
                {
                    cond_l = false;
                }
            }
        }
EndOfLoop:
    	return;
    }

    size_t tsp_ds::compute_idx(const size_t i, const size_t j, const size_t n) const
    {
        pagmo_assert( i!=j && i<n && j<n );
        return i*(n-1) + j - (j>i? 1:0);
    }

    void tsp_ds::compute_constraints_impl(constraint_vector &c, const decision_vector& x) const 
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
                    for (size_t j = 0; j < n_cities; j++) {
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
                for (size_t i = 0; i < n_cities; i++) {
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

    /// Computation of the phaseless distance between two planets
    /**
     * @parameter[in] i index of the first orbiting object
     * @parameter[in] j index of the second orbiting object
     * @return the distance between m_planets[i] and m_planets[j]
     */
    double tsp_ds::distance(decision_vector::size_type i, decision_vector::size_type j) const
    {
        // TODO: Edelbaum here instead? (see paper from Gatto-Casalino)
        using namespace std;
        kep_toolbox::array6D elements1 =  m_planets[i]->compute_elements(); /// Careful here what epoch? The default is 0 mjd2000, but the debris will be decayed by then and an exception may occur.
        kep_toolbox::array6D elements2 =  m_planets[j]->compute_elements();

        double a1 = elements1[0];
        double i1 = elements1[2];
        double W1 = elements1[3];
        double e1 = elements1[1];

        double a2 = elements2[0];
        double i2 = elements2[2];
        double W2 = elements2[3];
        double e2 = elements2[1];
        return three_impulses(a1,i1,W1,e1,a2,i2,W2,e2);
    }

    boost::array<int, 2> tsp_ds::compute_dimensions(decision_vector::size_type n_cities, base_tsp::encoding_type encoding)
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

    void tsp_ds::precompute_ephemerides() const 
    {
        kep_toolbox::array3D dumb0;
        std::vector<kep_toolbox::array3D> dumb1(m_planets.size(),dumb0);
        std::vector<std::vector<kep_toolbox::array3D> >dumb2(m_epochs.size(),dumb1);
        m_eph_r = dumb2;
        m_eph_v = dumb2;

        kep_toolbox::array6D dumb3;
        std::vector<kep_toolbox::array6D> dumb4(m_planets.size(),dumb3);
        std::vector<std::vector<kep_toolbox::array6D> >dumb5(m_epochs.size(),dumb4);
        m_eph_el = dumb5;
        // Precompute all ephemerides
        for (auto i = 0u; i < m_epochs.size(); ++i) 
        {
            for (auto j = 0u; j < m_planets.size(); ++j) {
                try
                {
                    m_planets[j]->eph(kep_toolbox::epoch(m_epochs[i], kep_toolbox::epoch::MJD2000), m_eph_r[i][j], m_eph_v[i][j]);
                    kep_toolbox::ic2par(m_eph_r[i][j], m_eph_v[i][j], m_mu, m_eph_el[i][j]);
                }
                catch( ... )
                {
                    std::cout << *m_planets[j] << std::endl;
                    std::cout << "At epoch: " << m_epochs[i] << std::endl;
                    std::cout << "planet idx: " << j << std::endl;
                    pagmo_throw(value_error, "Ephemerides computations caused an error (planet above)");
                }
            }
        }

    }

    void tsp_ds::compute_DVs(const decision_vector& tour, bool static_computations) const 
    {
        if (!static_computations)
        {
            for (auto i = 0u; i < get_n_cities() - 1;  ++i) 
            {
                m_DV[i] = distance_3_imp(tour[i], tour[(i + 1)], i);
            }
        } else {
            for (auto i = 0u; i < get_n_cities() - 1;  ++i) 
            {
                m_DV[i] = distance(tour[i], tour[(i + 1)]);
            }
        }
    }

    double tsp_ds::three_impulses(double a1, double i1, double W1, double e1, double a2, double i2, double W2, double e2) const
    {
        using namespace std;
        double Vi,Vf, DV1,DV2;

        // radius of apocenter starting orbit (km)
        double ra1 = a1 * (1. + e1);
        // radius of apocenter target orbit(km)
        double ra2 = a2 * (1. + e2);
        // relative inclination between orbits
        double cosiREL = cos(i1) * cos(i2) + sin(i1) * sin(i2) * cos(W1) * cos(W2) + sin(i1) * sin(i2) * sin(W1) * sin(W2);

        // radius of pericenter target orbit(km)
        double rp2 = a2 * (1. - e2);
        // radius of pericenter starting orbit (km)
        double rp1 = a1 * (1. - e1);

        if (ra1 > ra2) { // Strategy is Apocenter-Pericenter
            Vi = sqrt(m_mu * (2. / ra1 - 1. / a1));
            Vf = sqrt(m_mu * (2. / ra1 - 2. / (rp2 + ra1)));
            // Change Inclination + pericenter change
            DV1 = sqrt(Vi * Vi + Vf * Vf - 2. * Vi * Vf * cosiREL);
            // Apocenter Change
            DV2 = sqrt(m_mu) * abs(sqrt(2. / rp2 - 2. / (rp2 + ra1)) - sqrt(2. / rp2 - 2. / (rp2 + ra2)));
        }
        else {  // (ra1<ra2) Strategy is Pericenter-Apocenter
            // Apocenter Raise
            DV1 = sqrt(m_mu) * abs(sqrt(2. / rp1 - 2. / (rp1 + ra1)) - sqrt(2. / rp1 - 2. / (rp1 + ra2)));
            Vi = sqrt(m_mu * (2. / ra2 - 2. / (rp1 + ra2)));
            Vf = sqrt(m_mu * (2. / ra2 - 1. / a2));
            // Change Inclination + apocenter change
            DV2 = sqrt(abs(Vi * Vi + Vf * Vf - 2 * Vi * Vf * cosiREL));
        }
        return DV1 + DV2;
    }



    /// Computation of the distance between two planets using Lambert model
    /**
     * @parameter[in] dep index of the first orbiting object
     * @parameter[in] arr index of the second orbiting object
     * @parameter[in] ep_idx index of the epoch
     * @return the distance between m_planets[dep] and m_planets[arr]
     */

    double tsp_ds::distance_3_imp(const decision_vector::size_type dep, const decision_vector::size_type arr, const size_t ep_idx) const
    {
        // Computing the orbital parameters
        return three_impulses(m_eph_el[ep_idx][dep][0], m_eph_el[ep_idx][dep][2], m_eph_el[ep_idx][dep][3], m_eph_el[ep_idx][dep][1], m_eph_el[ep_idx][arr][0], m_eph_el[ep_idx][arr][2], m_eph_el[ep_idx][arr][3], m_eph_el[ep_idx][arr][1]);
    }


    /// Getter for m_planets
    /**
     * @return const reference to m_planets
     */
    const std::vector<kep_toolbox::planet::planet_ptr>& tsp_ds::get_planets() const
    { 
        return m_planets; 
    }

    /// Getter for m_max_DV
    /**
     * @return m_max_DV
     */
    double tsp_ds::get_max_DV() const
    { 
        return m_max_DV; 
    }


    /// Getter for m_values
    /**
     * @return const reference to m_values
     */
    const std::vector<double>&  tsp_ds::get_values() const
    { 
        return m_values; 
    }

    /// Getter for m_times
    /**
     * @return const reference to m_times
     */
    const decision_vector& tsp_ds::get_epochs() const
    { 
        return m_epochs; 
    }
    
    /// Getter for m_eph_el
    /**
     * @return const reference to m_eph_el
     */
    //const std::vector<std::vector<kep_toolbox::array6D>>& tsp_ds::get_eph_el() const
    //{ 
    //    return m_eph_el; 
    //}

    /// Returns the problem name
    std::string tsp_ds::get_name() const
    {
        return "Debris Selection TSP (TSP-DS)";
    }

    /// Extra human readable info for the problem.
    /**
     * @return a std::string containing a list of vertices and edges
     */
    std::string tsp_ds::human_readable_extra() const 
    {
        std::ostringstream oss;
        oss << "\n\tNumber of planets: " << get_n_cities() << '\n';
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
        oss << "\tObjects Names: " << std::endl << "\t\t[";
        for (auto i = 0u; i < m_planets.size(); ++i)
        {
            oss << m_planets[i]->get_name() << ", ";
            if (i == 4) break;
        }

        oss << "]" << std::endl << "\tEpochs: " << m_epochs << "\n";
        oss << "\ttObjects Values: " << m_values << std::endl;
        oss << "\tSpacecraft DV: " << m_max_DV << " [m/s]\n";
        
        return oss.str();
    }

    
}} //namespaces

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::tsp_ds)
