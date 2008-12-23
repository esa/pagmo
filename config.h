#ifndef PAGMO_CONFIG_H
#define PAGMO_CONFIG_H

#ifdef __GNUC__
#define GCC_VERSION (__GNUC__ * 100000 \
	+ __GNUC_MINOR__ * 1000 \
	+ __GNUC_PATCHLEVEL__ * 10)
#else
#error "The GCC compiler is needed to compile PaGMO."
#endif

#endif
