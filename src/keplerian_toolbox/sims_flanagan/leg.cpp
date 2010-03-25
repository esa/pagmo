#include "leg.h"
#include "sc_state.h"
#include "../astro_constants.h"
#include "../core_functions/array3D_operations.h"
#include "../core_functions/propagate_lagrangian.h"
#include"../../exceptions.h"
#include <vector>
#include <numeric>

namespace kep_toolbox{
    namespace sims_flanagan{


/// Overload the stream operator for kep_toolbox::sims_flanagan::leg
/**
 * Streams out the leg object in a human readable format
 *
 * \param[in] s stream to which the planet will be sent
 * \param[in] in leg to be sent to stream
 *
 * \return reference to s
 *
 */
	std::ostream &operator<<(std::ostream &s, const leg &in ){
	    s << "Number of segments: " << in.get_throttles_size() / 3 << std::endl;
	    s << "Departure date: " << in.get_t_i() << ", mjd2000: " << in.get_t_i().mjd2000() << std::endl;
	    s << "Arrival date: " << in.get_t_f() << ", mjd2000: " << in.get_t_f().mjd2000() << std::endl;
	    s << "Initial mass: " << in.get_x_i().get_mass() << " kg" << std::endl;
	    s << "Final mass: " << in.get_x_f().get_mass() << " kg" << std::endl;
	    s << "Absolute velocity at departure: " << norm(in.get_x_i().get_velocity()) << " m/s" << std::endl;
	    s << "Absolute velocity at arrival: " << norm(in.get_x_f().get_velocity()) << " m/s" << std::endl;

	    s << std::endl << "Throttles values: " << std::endl;
	    for (size_t i=0;i<in.get_throttles_size()/3;i++) {
		s << "\t\t\t" << in.throttles[i].get_value()[0] << " " << in.throttles[i].get_value()[1] << " " << in.throttles[i].get_value()[2] << std::endl;
	    }

	    std::vector<double> temp(in.get_throttles_size()/3);
	    in.get_throttles_con(temp.begin(), temp.end());
	    sc_state mism; in.get_mismatch_con(mism);
	    s << std::endl << "Mismatch at the midpoint: ";
	    std::cout << mism << std::endl;
	    s << "Throttle magnitude constraints: ";
	    for (size_t i=0;i< in.get_throttles_size()/3;i++) std::cout << temp[i] << " ";
	    std::cout << std::endl;
	    return s;
	}
    }
} //namespaces

/// Initialize a leg
/**
 * Initialize a leg assuming that the user has or will initialize separately the spacecraft
 * and the central body gravity parameter. The throttles are provided via two iterators pointing
 * to the beginning and to the end of a sequence of doubles containing the cartesian components
 * \f$ (x_1,y_1,z_1,x_2,y_2,z_2,...) \f$ of the throttles \f$ x_i,y_i,z_i \in [0,1]\f$. Needs to have dimension \f$ 3n \f$.
 *
 * \param[in] epoch_i Inital epoch for the leg
 * \param[in] state_i Initial sc_state (spacecraft state)
 * \param[in] throttles_start STL vector iterator pointing to the beginning of a cartesian throttle sequence. Throttles are numbers between 0 and 1.
 * \param[in] throttles_end STL vector iterator pointing to the end+1 of a cartesian throttle sequence. Throttles are numbers between 0 and 1.
 * \param[in] epoch_f Final epoch for the leg. Needs to be later than epoch_i
 * \param[in] state_f Final sc_state (spacecraft state)
 */

