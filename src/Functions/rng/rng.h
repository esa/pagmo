#ifndef PAGMO_RNG_H
#define PAGMO_RNG_H

#include <boost/random/lagged_fibonacci.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/thread/locks.hpp>
#include <boost/thread/mutex.hpp>
#include <ctime>

// This rng returns an unsigned integer in the [0,2**32-1] range.
typedef boost::mt19937 rng_uint32;
// This rng returns a double in the [0,1] range.
typedef boost::lagged_fibonacci607 rng_double;

// Generic thread-safe wrapper class around a Boost-like random number generator.
template <class Rng>
class static_rng {
	public:
		typedef typename Rng::result_type result_type;
		result_type operator()() {
			boost::lock_guard<boost::mutex> lock(m_mutex);
			return m_rng();
		}
	private:
		static boost::mutex	m_mutex;
		static Rng		m_rng;
};

template <class Rng>
boost::mutex static_rng<Rng>::m_mutex;

template <class Rng>
Rng static_rng<Rng>::m_rng(uint32_t(time(0)));

// Thread-safe uint32 rng.
typedef static_rng<rng_uint32> static_rng_uint32;

// Thread-safe double rng.
typedef static_rng<rng_double> static_rng_double;

#endif
