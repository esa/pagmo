#ifndef PAGMO_RNG_H
#define PAGMO_RNG_H

#include <boost/random/lagged_fibonacci.hpp>
#include <boost/random/mersenne_twister.hpp>

// This rng returns an unsigned integer in the [0,2**32-1] range.
typedef boost::mt19937 rng_uint32_type;
// This rng returns a double in the [0,1] range.
typedef boost::lagged_fibonacci607 rng_double_type;

#endif
