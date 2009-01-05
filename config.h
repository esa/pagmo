#ifndef PAGMO_CONFIG_H
#define PAGMO_CONFIG_H

#ifdef __GNUC__
#define GCC_VERSION (__GNUC__ * 100000 \
	+ __GNUC_MINOR__ * 1000 \
	+ __GNUC_PATCHLEVEL__ * 10)
#else
#error "The GCC compiler is needed to compile PaGMO."
#endif

#ifdef PAGMO_WIN32
#ifdef PAGMO_DLL_EXPORT_API
#define __PAGMO_VISIBLE __declspec(dllexport)
#elif defined ( PAGMO_DLL_IMPORT_API )
#define __PAGMO_VISIBLE __declspec(dllimport)
#else
#define __PAGMO_VISIBLE
#endif
#else
#define __PAGMO_VISIBLE __attribute__ ((visibility("default")))
#endif

#endif
