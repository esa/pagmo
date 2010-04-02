#include "planet.h"
#include"core_functions/ic2par.h"
#include"core_functions/par2ic.h"
#include"core_functions/convert_anomalies.h"
#include"exceptions.h"

namespace kep_toolbox{

    planet::planet(const epoch& ref_epoch, const array6D& orbital_elements_,
		   const double & mu_central_body_, const double & radius_, const double &mu_self_) :
	keplerian_elements(orbital_elements_),
	ref_mjd2000(ref_epoch.mjd2000()),
	radius(radius_), safe_radius(radius_*0.1),
	mu_self(mu_self_), mu_central_body(mu_central_body_),
	cached_epoch(std::numeric_limits<double>::quiet_NaN())
    {
	if (orbital_elements_.size() != 6) {
	    throw_value_error("The planet orbital element need to be a vector of dimension six");
	}
	if (orbital_elements_[0] <=0) {
	    throw_value_error("The planet semi-major axis needs to a positive number");
	}
	if (orbital_elements_[1] < 0 || orbital_elements_[1] >=1) {
	    throw_value_error("The planet eccentricity needs to be in [0,1)");
	}
	if (radius_ < 0) {
	    throw_value_error("The planet radius needs to be a positive number");
	}
	if (mu_central_body_ <= 0) {
	    throw_value_error("The central body gravitational parameter needs to be strictly positive");
	}
	if (mu_self_ < 0) {
	    throw_value_error("The gravitational parameter of the planet needs to be positive number");
	}
    }

    planet::planet(const std::string& name) :
    	cached_epoch(std::numeric_limits<double>::quiet_NaN())
    {
	if(name == "MERCURY")
	    initialize_planet(planet::MERCURY);
	else if(name == "VENUS")
	    initialize_planet(planet::VENUS);
	else if(name == "EARTH")
	    initialize_planet(planet::EARTH);
	else if(name == "MARS")
	    initialize_planet(planet::MARS);
	else if(name == "JUPITER")
	    initialize_planet(planet::JUPITER);
	else if(name == "SATURN")
	    initialize_planet(planet::SATURN);
	else if(name == "URANUS")
	    initialize_planet(planet::URANUS);
	else if(name == "NEPTUNE")
	    initialize_planet(planet::NEPTUNE);
	else
	    throw_value_error(std::string("unknown planet name ") + name);
    }
    
    planet::planet(const planet::common_name& name) :
    	cached_epoch(std::numeric_limits<double>::quiet_NaN())
    {
	initialize_planet(name);
    }
    
    void planet::initialize_planet(planet::common_name name) {
	ref_mjd2000 = 0;
	switch ( name ) {
	case planet::MERCURY: {
            double E[6] = {3.8717591e-01 * ASTRO_AU, 2.0563012e-01,7.0057713e+00 * ASTRO_DEG2RAD,4.8323635e+01 * ASTRO_DEG2RAD, 2.9131238e+01 * ASTRO_DEG2RAD, 1.7274971e+02 * ASTRO_DEG2RAD};
	    std::copy(E, E + 6, keplerian_elements.begin());
            radius = 2440000;
            safe_radius = radius + 200000;//radius*1.1;
            mu_self = 22032e9;
            mu_central_body = ASTRO_MU_SUN;
        }
	    break;
	case planet::VENUS: {
            double E[6] = {7.2347372e-01 * ASTRO_AU, 6.757291e-03, 3.3948519e+00 * ASTRO_DEG2RAD, 7.6659746e+01 * ASTRO_DEG2RAD, 5.5219953e+01 * ASTRO_DEG2RAD, 4.9298376e+01 * ASTRO_DEG2RAD};
	    std::copy(E, E + 6, keplerian_elements.begin());
            radius = 6052000;
            safe_radius = radius + 300000;//radius*1.1;
            mu_self = 324859e9;
            mu_central_body = ASTRO_MU_SUN;
        }
	    break;
	case planet::EARTH: {
            double E[6] = {9.9998805e-01 * ASTRO_AU, 1.6716812e-02, 8.8543531e-04 * ASTRO_DEG2RAD, 1.7540648e+02 * ASTRO_DEG2RAD, 2.8761578e+02 * ASTRO_DEG2RAD, 2.5760684e+02 * ASTRO_DEG2RAD};
	    std::copy(E, E + 6, keplerian_elements.begin());
            ref_mjd2000 = epoch::epoch(54000.0,epoch::MJD).mjd2000();
            radius = 6378000;
            safe_radius = radius*1.1;
            mu_self = 398600.4418e9;
            mu_central_body = ASTRO_MU_SUN;
        }
	    break;
	case planet::MARS: {
            double E[6] = { 1.5239844e+00 * ASTRO_AU, 9.3314935e-02, 1.8506136e+00 * ASTRO_DEG2RAD, 4.9535248e+01 * ASTRO_DEG2RAD, 2.865642e+02 * ASTRO_DEG2RAD, 1.909443e+01  * ASTRO_DEG2RAD };
	    std::copy(E, E + 6, keplerian_elements.begin());
            radius = 3397000;
            safe_radius = radius*1.1;
            mu_self = 42828e9;
            mu_central_body = ASTRO_MU_SUN;
        }
	    break;
	case planet::JUPITER: {
            double E[6] = {5.2107645e+00 * ASTRO_AU, 4.9715759e-02, 1.3044197e+00 * ASTRO_DEG2RAD, 1.0044249e+02 * ASTRO_DEG2RAD,2.7550414e+02 * ASTRO_DEG2RAD,1.8387562e+01 * ASTRO_DEG2RAD};
	    std::copy(E, E + 6, keplerian_elements.begin());
            radius = 71492000;
            safe_radius = 671492000;
            mu_self = 126686534e9;
            mu_central_body = ASTRO_MU_SUN;
        }
	    break;
	case planet::SATURN: {
            double E[6] = {9.5869202e+00 * ASTRO_AU, 5.5944004e-02, 2.4847848e+00 * ASTRO_DEG2RAD, 1.1361884e+02 * ASTRO_DEG2RAD, 3.3583259e+02 * ASTRO_DEG2RAD,-3.9463561e+01 * ASTRO_DEG2RAD};
	    std::copy(E, E + 6, keplerian_elements.begin());
            radius = 60330000;
            safe_radius = radius*1.1;
            mu_self = 37931187e9;
            mu_central_body = ASTRO_MU_SUN;
        }
	    break;
	case planet::URANUS: {
            double E[6] = {1.9234105e+01 * ASTRO_AU, 4.4369076e-02, 7.7287008e-01 * ASTRO_DEG2RAD, 7.3908932e+01 * ASTRO_DEG2RAD, 9.6656163e+01 * ASTRO_DEG2RAD, 1.4291587e+02 * ASTRO_DEG2RAD};
	    std::copy(E, E + 6, keplerian_elements.begin());
            radius = 25559000;
            safe_radius = radius*1.1;
            mu_self = 5793939e9;
            mu_central_body = ASTRO_MU_SUN;
        }
	    break;
	case planet::NEPTUNE: {
            double E[6] = {3.0111359e+01 * ASTRO_AU, 1.1211871e-02, 1.7672166e+00 * ASTRO_DEG2RAD, 1.3176686e+02 * ASTRO_DEG2RAD, 2.6539295e+02 * ASTRO_DEG2RAD, -9.1954294e+01 * ASTRO_DEG2RAD};
	    std::copy(E, E + 6, keplerian_elements.begin());
            radius = 24764000;
            safe_radius = radius*1.1;
            mu_self = 6836528e9;
            mu_central_body = ASTRO_MU_SUN;
        }
	    break;
	}
	time_coefficient = sqrt(mu_central_body / pow(keplerian_elements[0],3));
    }
    
    void planet::get_eph(const epoch& when, array3D &r, array3D &v) const{
	if(when.mjd2000() != cached_epoch.mjd2000()) {
	    double elements[6];
	    std::copy(keplerian_elements.begin(), keplerian_elements.end(), elements);
	    double dt = (when.mjd2000() - ref_mjd2000) * ASTRO_DAY2SEC;
	    elements[5] += time_coefficient * dt;
	    elements[5] = m2e(elements[5],elements[1]);
	    par2ic(elements, mu_central_body, cached_r, cached_v);
	    cached_epoch = when;
	}
	
	r = cached_r;
	v = cached_v;
    }

    array3D planet::get_position(const epoch& when) const {
	array3D r, v;
	get_eph(when, r, v);
	return (r);
    }

    array3D planet::get_velocity(const epoch& when) const{
	array3D r, v;
	get_eph(when, r, v);
	return (v);
    }
    
    array6D planet::get_elements(const epoch& when) const{
	array6D elements(keplerian_elements);
	double dt = (when.mjd2000() - ref_mjd2000) * ASTRO_DAY2SEC;
	elements[5] += time_coefficient * dt;
	return ( elements );
    }

}
std::ostream &kep_toolbox::operator<<(std::ostream &s, const kep_toolbox::planet &body) {
    s << "Keplerian Planet: " << std::endl;
    s << "Own gravity parameter: " << body.mu_self << std::endl;
    s << "Central body gravity parameter: " << body.mu_central_body << std::endl;
    s << "Planet radius: " << body.radius << std::endl;
    s << "Planet keplerian elements: "<<std::endl;
    array6D elem = body.get_elements(epoch(body.ref_mjd2000));
    s << "Semi major axis (AU): " << elem[0] / ASTRO_AU << std::endl;
    s << "Eccentricity: " << elem[1] << std::endl;
    s << "Inclination (deg.): " << elem[2] * ASTRO_RAD2DEG << std::endl;
    s << "Big Omega (deg.): " << elem[3] * ASTRO_RAD2DEG << std::endl;
    s << "Small omega (deg.): " << elem[4] * ASTRO_RAD2DEG << std::endl;
    s << "Mean anomaly (deg.): " << elem[5] * ASTRO_RAD2DEG << std::endl;
    s << "Elements reference epoch: " << epoch(body.ref_mjd2000) << std::endl;
    return s;
}
