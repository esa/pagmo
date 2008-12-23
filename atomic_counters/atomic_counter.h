#ifndef PAGMO_ATOMIC_COUNTER_H
#define PAGMO_ATOMIC_COUNTER_H

#if defined( __GNUC__ ) && GCC_VERSION >= 401000

#include "atomic_counter_gcc_41.h"

namespace PaGMO
{
	typedef atomic_counter_gcc_41<size_t> atomic_counter_size_t;
}

#else // Not GCC or GCC < 4.1.

// TODO: for MSVC, use its atomic builtins instead of the generic counter.
#include "atomic_counter_generic.h"

namespace PaGMO
{
	typedef atomic_counter_generic<size_t> atomic_counter_size_t;
}

#endif // Compiler selection in case of MT.

#endif
