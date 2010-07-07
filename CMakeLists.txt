PROJECT(PaGMO)

ENABLE_TESTING()

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.0)
# Setup and include CMake modules.
SET(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake_modules)

INCLUDE(CheckTypeSize)
INCLUDE(PaGMOCPackSetup)

# Initialise (empty) list of libraries to link against.
SET(MANDATORY_LIBRARIES "")

# Configuration for GCC.
IF(CMAKE_COMPILER_IS_GNUCXX)
	# Add to the base flags extra warnings, and set up message reporting from the compiler so that it is Eclipse-friendly.
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wdisabled-optimization -fmessage-length=0 -Wnon-virtual-dtor")
	# Additional flags to turn off some GCC warnings that in practice clutter the compilation output.
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-missing-field-initializers -fno-strict-aliasing -Wno-deprecated")
	# Suggested for multithreaded code.
	ADD_DEFINITIONS(-D_REENTRANT)
	# Atomic builtins for GCC.
	TRY_COMPILE(WORKING_ATOMIC_BUILTINS ${CMAKE_BINARY_DIR}/compile_tests/ ${CMAKE_SOURCE_DIR}/cmake_modules/gcc_atomic_builtins_test.cpp COMPILE_DEFINITIONS ${CMAKE_CCXX_FLAGS})
	IF(WORKING_ATOMIC_BUILTINS)
		MESSAGE(STATUS "GCC is correctly setup to support atomic builtins.")
		ADD_DEFINITIONS(-D_PAGMO_GCC_ATOMIC_BUILTINS)
	ELSE(WORKING_ATOMIC_BUILTINS)
		MESSAGE(STATUS "Either this GCC version does not support atomic builtins or the CXXFLAGS are not properly set to enable them.")
		MESSAGE(STATUS "Please note that atomic builtins are supported from version 4.1.0 of GCC and need an appropriate '-march' flag (e.g., at least '-march=i486' on x86 architectures).")
	ENDIF(WORKING_ATOMIC_BUILTINS)
ENDIF(CMAKE_COMPILER_IS_GNUCXX)

# Configuration for MSVC.
IF(MSVC)
	MESSAGE(STATUS "Enabling MSVC atomic builtins.")
	ADD_DEFINITIONS(-D_PAGMO_MSVC_ATOMIC_BUILTINS)
ENDIF(MSVC)

# Set default build type to "Release", change it in the GUI if you need to build with debug.
IF(NOT CMAKE_BUILD_TYPE)
	SET(CMAKE_BUILD_TYPE Release CACHE STRING
		"Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
	FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

# Use CMake's 2.6 new policy for library paths.
IF(COMMAND CMAKE_POLICY)
    CMAKE_POLICY(SET CMP0003 NEW)
ENDIF(COMMAND CMAKE_POLICY)

# Provides build options to CMake.
# Build Option: when active the file main.cpp is built and linked to the PaGMO static library.
OPTION(BUILD_MAIN "Build 'main.cpp'." ON)

# Build Option: when active the bindings to Python are compiled and linked (installation).
OPTION(BUILD_PYGMO "Build Python bindings." OFF)

# Build Option: minimisers from the GNU scientific library (GSL).
OPTION(ENABLE_GSL "Enable support for GSL minimisers (requires GSL >= 1.13)." OFF)

# Build Option: algorithms from the NLopt library.
OPTION(ENABLE_NLOPT "Enable support for NLopt minimisers." OFF)

# Build Option: SNOPT optimiser.
OPTION(ENABLE_SNOPT "Enable support for SNOPT minimiser." OFF)

# Build Option: IPOPT optimiser.
OPTION(ENABLE_IPOPT "Enable support for IPOPT minimiser." OFF)

# Build option: enable test set.
OPTION(ENABLE_TESTS "Build test set." ON)

# Build Option: build executable for the examples
OPTION(BUILD_EXAMPLES "Creates executables for all scripts in folder examples" OFF)

# If GSL support is requested, look for the library.
IF(ENABLE_GSL)
	FIND_LIBRARY(GSL_GSL_LIBRARY NAMES gsl)
	FIND_LIBRARY(GSL_GSLCBLAS_LIBRARY NAMES gslcblas)
	FIND_PATH(GSL_INCLUDE_DIR NAMES gsl/gsl_version.h)
	MESSAGE(STATUS "GSL library: ${GSL_GSL_LIBRARY}")
	MESSAGE(STATUS "GSL CBLAS library: ${GSL_GSLCBLAS_LIBRARY}")
	MESSAGE(STATUS "GSL include dir: ${GSL_INCLUDE_DIR}")
	IF(NOT GSL_GSL_LIBRARY OR NOT GSL_INCLUDE_DIR OR NOT GSL_GSLCBLAS_LIBRARY)
		MESSAGE(FATAL_ERROR "GSL support was requested, but at least one of the required libraries (gsl and gslcblas) could not be located.")
	ENDIF(NOT GSL_GSL_LIBRARY OR NOT GSL_INCLUDE_DIR OR NOT GSL_GSLCBLAS_LIBRARY)
	MESSAGE(STATUS "Running the GSL version checker.")
	TRY_RUN(SUFFICIENT_GSL_VERSION_RUN SUFFICIENT_GSL_VERSION_COMPILE ${CMAKE_BINARY_DIR}/compile_tests/ ${CMAKE_SOURCE_DIR}/cmake_modules/gsl_version.cpp
		COMPILE_DEFINITIONS "${CMAKE_CCXX_FLAGS} -I\"${CMAKE_SOURCE_DIR}/external_headers\" -I\"${GSL_INCLUDE_DIR}\"" RUN_OUTPUT_VARIABLE SUFFICIENT_GSL_VERSION_OUTPUT)
	MESSAGE(STATUS "GSL version checker output: ${SUFFICIENT_GSL_VERSION_OUTPUT}")
	IF(NOT SUFFICIENT_GSL_VERSION_COMPILE)
		MESSAGE(FATAL_ERROR "Error compiling the GSL version checker.")
	ENDIF(NOT SUFFICIENT_GSL_VERSION_COMPILE)
	IF(NOT ${SUFFICIENT_GSL_VERSION_RUN} EQUAL 0)
		MESSAGE(FATAL_ERROR "Error running the GSL version checker: either the GSL version is < 1.13 or your GSL installation is broken.")
	ELSE(NOT ${SUFFICIENT_GSL_VERSION_RUN} EQUAL 0)
		MESSAGE(STATUS "GSL version is fine.")
	ENDIF(NOT ${SUFFICIENT_GSL_VERSION_RUN} EQUAL 0)
	SET(MANDATORY_LIBRARIES ${MANDATORY_LIBRARIES} ${GSL_GSL_LIBRARY} ${GSL_GSLCBLAS_LIBRARY})
	ADD_DEFINITIONS(-DPAGMO_ENABLE_GSL)
	INCLUDE_DIRECTORIES("${GSL_INCLUDE_DIR}")
ENDIF(ENABLE_GSL)

# If NLopt support is requested, look for the library.
IF(ENABLE_NLOPT)
	FIND_LIBRARY(NLOPT_LIBRARY NAMES nlopt_cxx)
	FIND_PATH(NLOPT_INCLUDE_DIR NAMES nlopt.h)
	IF(NOT NLOPT_LIBRARY OR NOT NLOPT_INCLUDE_DIR)
		MESSAGE(FATAL_ERROR "NLopt support was requested, but the library could not be located. Please make sure that NLopt was compiled with C++ support (the library will be called 'nlopt_cxx')")
	ENDIF(NOT NLOPT_LIBRARY OR NOT NLOPT_INCLUDE_DIR)
	MESSAGE(STATUS "NLopt library: ${NLOPT_LIBRARY}")
	MESSAGE(STATUS "NLopt include dir: ${NLOPT_INCLUDE_DIR}")
	SET(MANDATORY_LIBRARIES ${MANDATORY_LIBRARIES} ${NLOPT_LIBRARY})
	ADD_DEFINITIONS(-DPAGMO_ENABLE_NLOPT)
	INCLUDE_DIRECTORIES("${NLOPT_INCLUDE_DIR}")
ENDIF(ENABLE_NLOPT)

# If SNOPT support is requested, look for the libraries.
IF(ENABLE_SNOPT)
	FIND_LIBRARY(SNOPT_SNOPT_LIBRARY NAMES snopt)
	IF(NOT SNOPT_SNOPT_LIBRARY)
		MESSAGE(FATAL_ERROR "SNOPT support was requested, but the snopt library could not be located.")
	ENDIF(NOT SNOPT_SNOPT_LIBRARY)
	MESSAGE(STATUS "snopt library: ${SNOPT_SNOPT_LIBRARY}")

	FIND_LIBRARY(SNOPT_SNPRINT_LIBRARY NAMES snprint)
	IF(NOT SNOPT_SNPRINT_LIBRARY)
		MESSAGE(FATAL_ERROR "SNOPT support was requested, but the snprint library could not be located.")
	ENDIF(NOT SNOPT_SNPRINT_LIBRARY)
	MESSAGE(STATUS "snprint library: ${SNOPT_SNPRINT_LIBRARY}")

	FIND_LIBRARY(SNOPT_BLAS_LIBRARY NAMES blas)
	IF(NOT SNOPT_BLAS_LIBRARY)
		MESSAGE(FATAL_ERROR "SNOPT support was requested, but the blas library could not be located.")
	ENDIF(NOT SNOPT_BLAS_LIBRARY)
	MESSAGE(STATUS "blas library: ${SNOPT_BLAS_LIBRARY}")

	FIND_LIBRARY(SNOPT_F2C_LIBRARY NAMES f2c)
	IF(NOT SNOPT_F2C_LIBRARY)
		MESSAGE(FATAL_ERROR "SNOPT support was requested, but the f2c library could not be located.")
	ENDIF(NOT SNOPT_F2C_LIBRARY)
	MESSAGE(STATUS "f2c library: ${SNOPT_F2C_LIBRARY}")

	FIND_LIBRARY(SNOPT_M_LIBRARY NAMES m)
	IF(NOT SNOPT_M_LIBRARY)
		MESSAGE(FATAL_ERROR "SNOPT support was requested, but the m library could not be located.")
	ENDIF(NOT SNOPT_M_LIBRARY)
	MESSAGE(STATUS "m library: ${SNOPT_M_LIBRARY}")

	SET(MANDATORY_LIBRARIES ${MANDATORY_LIBRARIES} ${SNOPT_SNOPT_LIBRARY} ${SNOPT_SNPRINT_LIBRARY} ${SNOPT_BLAS_LIBRARY} ${SNOPT_F2C_LIBRARY} ${SNOPT_M_LIBRARY})

	#If gfortran is found the library is included as last in case snopt has been compiled using gfortran
	FIND_LIBRARY(SNOPT_GFORTRAN_LIBRARY NAMES gfortran)
	IF(SNOPT_GFORTRAN_LIBRARY)
		SET(MANDATORY_LIBRARIES ${MANDATORY_LIBRARIES} ${SNOPT_GFORTRAN_LIBRARY})
		MESSAGE(STATUS "gfortran library was found, linking in.")
	ELSE(SNOPT_GFORTRAN_LIBRARY)
		MESSAGE(STATUS "gfortran library was not found.")
	ENDIF(SNOPT_GFORTRAN_LIBRARY)
	ADD_DEFINITIONS(-DPAGMO_ENABLE_SNOPT)
ENDIF(ENABLE_SNOPT)

# If IPOPT support is requested, look for the libraries.
IF(ENABLE_IPOPT)
	FIND_LIBRARY(IPOPT_LIBRARY NAMES ipopt)
	FIND_LIBRARY(LAPACK_LIBRARY NAMES lapack)
	FIND_PATH(IPOPT_INCLUDE_DIR NAMES coin/IpTNLP.hpp)
	IF(NOT IPOPT_LIBRARY OR NOT IPOPT_INCLUDE_DIR OR NOT LAPACK_LIBRARY)
		MESSAGE(FATAL_ERROR "IPOPT support was requested, but (one of) the libraries could not be located. Please make sure that you have ipopt library and the headers installed in the default coin/ directory")
	ENDIF(NOT IPOPT_LIBRARY OR NOT IPOPT_INCLUDE_DIR OR NOT LAPACK_LIBRARY)
	MESSAGE(STATUS "IPOPT library: ${IPOPT_LIBRARY}")
	MESSAGE(STATUS "IPOPT include dir: ${IPOPT_INCLUDE_DIR}")
	SET(MANDATORY_LIBRARIES ${MANDATORY_LIBRARIES} ${IPOPT_LIBRARY} ${LAPACK_LIBRARY})
	ADD_DEFINITIONS(-DPAGMO_ENABLE_IPOPT)
	INCLUDE_DIRECTORIES("${IPOPT_INCLUDE_DIR}")
ENDIF(ENABLE_IPOPT)

# Let's include Boost's headers.
INCLUDE_DIRECTORIES("${CMAKE_SOURCE_DIR}/external_headers")

# Some common platform switches.
SET(PAGMO_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX})
IF(UNIX)
	# Install path for libraries.
	SET(LIB_INSTALL_PATH "lib")
	# NOTICE: use FindThreads library from CMake here instead?
	# If we are in Unix, let's require pthread.
	FIND_LIBRARY(PTHREAD_LIBRARY pthread)
	IF(NOT PTHREAD_LIBRARY)
		MESSAGE(FATAL_ERROR "POSIX thread library not found, please install it or specify manually where it is located.")
	ENDIF(NOT PTHREAD_LIBRARY)
	SET(MANDATORY_LIBRARIES ${MANDATORY_LIBRARIES} pthread)
	# Enable the pthread flag in Unix only if the compiler is GNU.
	IF(CMAKE_COMPILER_IS_GNUCXX)
		SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
	ENDIF(CMAKE_COMPILER_IS_GNUCXX)
ELSE(UNIX)
	IF(NOT WIN32)
		MESSAGE(FATAL_ERROR "Sorry, at the moment pagmo can be built only on Unix and Win32 environments.")
	ENDIF(NOT WIN32)
	SET(LIB_INSTALL_PATH .)
	SET(PAGMO_INSTALL_PREFIX "${PAGMO_INSTALL_PREFIX} ${PAGMO_VERSION}")
	ADD_DEFINITIONS(-DPAGMO_WIN32)
	IF(MINGW)
		# Linker and compiler flags needed for the support of multithread applications.
		SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS_INIT} -mthreads")
		SET(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS_INIT} -mthreads")
		SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS_INIT} -mthreads")
		SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mthreads")
	ENDIF(MINGW)
ENDIF(UNIX)

CHECK_TYPE_SIZE("void *" POINTER_SIZE)
MESSAGE(STATUS "Pointer size = ${POINTER_SIZE}")

IF(BUILD_PYGMO)
	INCLUDE(PaGMOPythonSetup)
ENDIF(BUILD_PYGMO)

ADD_SUBDIRECTORY("${CMAKE_SOURCE_DIR}/src")
IF(BUILD_PYGMO)
	ADD_SUBDIRECTORY("${CMAKE_SOURCE_DIR}/PyGMO")
	IF(WIN32)
		SET(WIN_PKG_ROOT ${CMAKE_SOURCE_DIR}/windows_package)
		# Copy license file into windows package directory and install it.
		CONFIGURE_FILE("${CMAKE_SOURCE_DIR}/COPYING" "${WIN_PKG_ROOT}/license.txt")
		CONFIGURE_FILE("${CMAKE_SOURCE_DIR}/README" "${WIN_PKG_ROOT}/readme.txt")
		INSTALL(FILES "${WIN_PKG_ROOT}/license.txt" "${WIN_PKG_ROOT}/readme.txt" DESTINATION .)
		# Install Console.
		INSTALL(FILES
			${WIN_PKG_ROOT}/console/Console.exe
			${WIN_PKG_ROOT}/console/FreeImage.dll
			${WIN_PKG_ROOT}/console/console.chm
			${WIN_PKG_ROOT}/console/console.xml
			${WIN_PKG_ROOT}/console/ConsoleHook.dll
			${WIN_PKG_ROOT}/console/FreeImagePlus.dll
			DESTINATION .
		)
		INSTALL(FILES
			${WIN_PKG_ROOT}/console/Microsoft.VC90.CRT/msvcm90.dll
			${WIN_PKG_ROOT}/console/Microsoft.VC90.CRT/msvcr90.dll
			${WIN_PKG_ROOT}/console/Microsoft.VC90.CRT/msvcp90.dll
			${WIN_PKG_ROOT}/console/Microsoft.VC90.CRT/Microsoft.VC90.CRT.manifest
			DESTINATION Microsoft.VC90.CRT
		)
		# Install IPython profile file.
		INSTALL(FILES
			${CMAKE_SOURCE_DIR}/extras/ipy_profile_pygmo.py
			DESTINATION .
		)
	ENDIF(WIN32)
ENDIF(BUILD_PYGMO)

# Link main to pagmo_static library.
IF(BUILD_MAIN)
	ADD_EXECUTABLE(main main.cpp)
	TARGET_LINK_LIBRARIES(main  ${MANDATORY_LIBRARIES} pagmo_static)
ENDIF(BUILD_MAIN)

IF(ENABLE_TESTS)
	ADD_SUBDIRECTORY("${CMAKE_SOURCE_DIR}/tests")
ENDIF(ENABLE_TESTS)

IF(BUILD_EXAMPLES)
	ADD_SUBDIRECTORY("${CMAKE_SOURCE_DIR}/examples")
ENDIF(BUILD_EXAMPLES)

INCLUDE(CPack)
