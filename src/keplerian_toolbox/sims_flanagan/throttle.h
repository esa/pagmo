#ifndef THROTTLE_H
#define THROTTLE_H

namespace kep_toolbox {
/// Objects dealing with the Sims-Flanagan low-thrus trajectory model
/**
 * This namespace contains the routines that allow building and evaluating low-thrust trajectories using the
 * Sims-Flanagan type of approach.
 */
    namespace sims_flanagan{

	class throttle {
	public:
	    throttle() {}
	    
	    throttle(epoch _start, epoch _end, const array3D& _value)
		: start(_start), end(_end), value(_value) {}

	    epoch get_start() const {
		return start;
	    }
	    
	    epoch get_end() const {
		return end;
	    }

	    const array3D& get_value() const {
		return value;
	    }
	    
	private:
	    epoch start;
	    epoch end;
	    array3D value;
	};
    };
};

#endif
