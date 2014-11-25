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

#include "tsp_ads.h"
#include "../exceptions.h"
#include "../population.h"
#include "../keplerian_toolbox/astro_constants.h"
#include "../keplerian_toolbox/lambert_problem.h"
#include "../keplerian_toolbox/core_functions/array3D_operations.h"
#include "../keplerian_toolbox/epoch.h"


namespace pagmo { namespace problem {

       /// Constructor
    /**
     * Constructs a TSP - Asteroids-Debris Selection Problem from 
     * path length and the chosen encoding
     *
     * @param[in] planets         an std::vector of kep_toolbox::planet_ptr containing the orbiting objects one wants to visit
     * @param[in] values          a std::vector representing the value (score) of a visit (is the same order as planets)
     * @param[in] max_DV          the maximum DV available on board the spacecraft (corresponding to max_path_length parameter in the problem::tsp_cs problem)
     * @param[in] times           an std::vector of doubles containing [t0,T1,T2,...] the starting mission epoch and the time of flights between orbiting objects
     * @param[in] waiting_time    waiting time at one orbiting object before the next transfer starts (i.e. t_0+T_1+waiting_time+.....) 
     * @param[in] encoding        a pagmo::problem::tsp::encoding representing the chosen encoding
     */
    tsp_ads::tsp_ads (            
        const std::vector<kep_toolbox::planet_ptr>& planets, 
        const std::vector<double>& values,
        const double max_DV, 
        const std::vector<double>&  epochs, 
        const double waiting_time, 
        const base_tsp::encoding_type & encoding
        ) : base_tsp(
                planets.size(), 
                compute_dimensions(planets.size(), encoding)[0],
                compute_dimensions(planets.size(), encoding)[1],
                encoding
        ),  m_planets(planets), m_values(values), m_max_DV(max_DV), m_epochs(epochs), m_waiting_time(waiting_time)
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

        //TODO: check that all planets have the same mu

        m_min_value = *std::min(m_values.begin(), m_values.end());
    }

    /// Clone method.
    base_ptr tsp_ads::clone() const
    {
        return base_ptr(new tsp_ads(*this));
    }

    boost::array<int, 2> tsp_ads::compute_dimensions(decision_vector::size_type n_cities, base_tsp::encoding_type encoding)
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

    void tsp_ads::objfun_impl(fitness_vector &f, const decision_vector& x) const 
    {
        f[0]=0;
        decision_vector tour;
        decision_vector::size_type n_cities = get_n_cities();
        double cum_p, saved_length;
        decision_vector::size_type dumb1, dumb2;

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
        find_best_selection(tour, cum_p, saved_length, dumb1, dumb2);
        f[0] = -(cum_p + (1 - m_min_value) * n_cities + saved_length / m_max_DV);
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
     * @param[out] retval_l the total cumulative length of the subpath
     * @param[out] retval_it_l the id of the city where the subpath starts
     * @param[out] retval_it_r the id of the city where the subpath ends
     * @throws value_error if the input tour length is not equal to the city number

     */
    void tsp_ads::find_best_selection(const decision_vector& tour, double& retval_p, double& retval_l, decision_vector::size_type& retval_it_l, decision_vector::size_type& retval_it_r) const
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

        // Main body of the double loop
        while(cond_l)
        {
            while(cond_r) 
            {
                // We increment the right "pointer" updating the value and length of the path
                saved_length -= distance_lambert(tour[it_r % n_cities], tour[(it_r + 1) % n_cities]);
                cum_p += m_values[(it_r + 1) % n_cities];
                it_r += 1;

                // We update the various retvals only if the new subpath is valid
                if (saved_length < 0 || (it_l % n_cities == it_r % n_cities))
                {
                    cond_r = false;
                }
                else if (cum_p > retval_p)
                {
                    retval_p = cum_p;
                    retval_l = saved_length;
                    retval_it_l = it_l % n_cities;
                    retval_it_r = it_r % n_cities;
                }
                else if (cum_p == retval_p)
                {
                    if (saved_length > retval_l)
                    {
                        retval_p = cum_p;
                        retval_l = saved_length;
                        retval_it_l = it_l % n_cities;
                        retval_it_r = it_r % n_cities;
                    }
                }

                if (it_r==n_cities-1) 
                {
                    goto EndOfLoop;
                }
            }
            // We get out if all cities are included in the current path
            if (it_l % n_cities == it_r % n_cities)
            {
                cond_l = false;
            }
            else
            {
                // We increment the left "pointer" updating the value and length of the path
                saved_length += distance_lambert(tour[it_l % n_cities], tour[(it_l + 1) % n_cities]);
                cum_p -= m_values[it_l];
                it_l += 1;
                // We update the various retvals only if the new subpath is valid
                if (saved_length > 0)
                {
                    cond_r = true;
                    if (cum_p > retval_p)
                    {
                        retval_p = cum_p;
                        retval_l = saved_length;
                        retval_it_l = it_l % n_cities;
                        retval_it_r = it_r % n_cities;
                    }
                    else if (cum_p == retval_p)
                    {
                        if (saved_length > retval_l)
                        {
                            retval_p = cum_p;
                            retval_l = saved_length;
                            retval_it_l = it_l % n_cities;
                            retval_it_r = it_r % n_cities;
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

    size_t tsp_ads::compute_idx(const size_t i, const size_t j, const size_t n) const
    {
        pagmo_assert( i!=j && i<n && j<n );
        return i*(n-1) + j - (j>i? 1:0);
    }

    void tsp_ads::compute_constraints_impl(constraint_vector &c, const decision_vector& x) const 
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
    double tsp_ads::distance(decision_vector::size_type i, decision_vector::size_type j) const
    {
        // TODO: Edelbaum here instead? (see paper from Gatto-Casalino)
        using namespace std;
        double DV1=0, DV2=0, Vi=0, Vf=0;
        double MU = m_planets[i]->get_mu_central_body();
        kep_toolbox::array6D elements1 =  m_planets[i]->get_elements();
        kep_toolbox::array6D elements2 =  m_planets[j]->get_elements();

        double a1 = elements1[0];
        double i1 = elements1[2];
        double W1 = elements1[3];
        double e1 = elements1[1];

        double a2 = elements2[0];
        double i2 = elements2[2];
        double W2 = elements2[3];
        double e2 = elements2[1];
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
            Vi = sqrt(MU * (2. / ra1 - 1. / a1));
            Vf = sqrt(MU * (2. / ra1 - 2. / (rp2 + ra1)));
            // Change Inclination + pericenter change
            DV1 = sqrt(Vi * Vi + Vf * Vf - 2. * Vi * Vf * cosiREL);
            // Apocenter Change
            DV2 = sqrt(MU) * abs(sqrt(2. / rp2 - 2. / (rp2 + ra1)) - sqrt(2. / rp2 - 2. / (rp2 + ra2)));
        }
        else {  // (ra1<ra2) Strategy is Pericenter-Apocenter
            // Apocenter Raise
            DV1 = sqrt(MU) * abs(sqrt(2. / rp1 - 2. / (rp1 + ra1)) - sqrt(2. / rp1 - 2. / (rp1 + ra2)));
            Vi = sqrt(MU * (2. / ra2 - 2. / (rp1 + ra2)));
            Vf = sqrt(MU * (2. / ra2 - 1. / a2));
            // Change Inclination + apocenter change
            DV2 = sqrt(abs(Vi * Vi + Vf * Vf - 2 * Vi * Vf * cosiREL));
        }
        return DV1 + DV2;
    }

    /// Computation of the distance between two planets using Lambert model
    /**
     * @parameter[in] i index of the first orbiting object
     * @parameter[in] j index of the second orbiting object
     * @return the distance between m_planets[i] and m_planets[j]
     */
    double tsp_ads::distance_lambert(decision_vector::size_type dep, decision_vector::size_type arr) const
    {
        using namespace std;
        double t1 = m_epochs[dep];
        double t2 = m_epochs[arr];
        if (t2<=t1) 
        {
            pagmo_throw(value_error, "Epoch at arrival smaller than epoch at departure");
        }
        kep_toolbox::array3D r1,v1,r2,v2,DV1,DV2;

        //Ephemerides computations
        m_planets[dep]->get_eph(kep_toolbox::epoch(t1, kep_toolbox::epoch::MJD2000), r1, v1);
        m_planets[arr]->get_eph(kep_toolbox::epoch(t2, kep_toolbox::epoch::MJD2000), r2, v2);

        // Lambert's problem solution
        kep_toolbox::lambert_problem l(r1, r1, (t2 - t1) * ASTRO_DAY2SEC, m_planets[dep]->get_mu_central_body(), 0, 0);
        kep_toolbox::diff(DV1, l.get_v1()[0], v1);
        kep_toolbox::diff(DV2, l.get_v2()[0], v2);
        
        return kep_toolbox::norm(DV1) + kep_toolbox::norm(DV2);
    }

    /// Getter for m_planets
    /**
     * @return const reference to m_planets
     */
    const std::vector<kep_toolbox::planet_ptr>& tsp_ads::get_planets() const
    { 
        return m_planets; 
    }

    /// Getter for m_max_DV
    /**
     * @return m_max_DV
     */
    double tsp_ads::get_max_DV() const
    { 
        return m_max_DV; 
    }


    /// Getter for m_values
    /**
     * @return const reference to m_values
     */
    const std::vector<double>&  tsp_ads::get_values() const
    { 
        return m_values; 
    }

    /// Getter for m_times
    /**
     * @return const reference to m_times
     */
    const decision_vector& tsp_ads::get_epochs() const
    { 
        return m_epochs; 
    }

    /// Getter for m_waiting_time
    /**
     * @return m_waiting_time
     */
    double tsp_ads::get_waiting_time() const
    { 
        return m_waiting_time; 
    }

    /// Returns the problem name
    std::string tsp_ads::get_name() const
    {
        return "Asteroids / Debris Selection TSP (TSP-ADS)";
    }

    /// Extra human readable info for the problem.
    /**
     * @return a std::string containing a list of vertices and edges
     */
    std::string tsp_ads::human_readable_extra() const 
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
        oss << "\tPlanets Names: " << std::endl << "\t\t[";
        for (auto i = 0u; i < m_planets.size(); ++i)
        {
            oss << m_planets[i]->get_name() << ", ";
            if (i == 4) break;
        }

        oss << "]" << std::endl << "\tEpochs: " << m_epochs << "\n";
        oss << "\tPlanets Values: " << m_values << std::endl;
        oss << "\tSpacecraft DV: " << m_max_DV << " [m/s]\n";
        

        return oss.str();
    }

    
}} //namespaces

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::tsp_ads)