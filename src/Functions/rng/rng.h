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

template <class Rng>
class mt_rng {
	public:
		typedef typename Rng::result_type result_type;
		mt_rng():m_mutex(),m_rng(uint32_t(time(0))) {}
		result_type operator()() {
			boost::lock_guard<boost::mutex> lock(m_mutex);
			return m_rng();
		}
	private:
		boost::mutex	m_mutex;
		Rng		m_rng;
};

typedef mt_rng<rng_double> mt_rng_double;

#endif
