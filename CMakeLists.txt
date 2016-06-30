PROJECT(PaGMO)

# We enable cmake test capability
ENABLE_TESTING()

# We set the minimum required cmake version
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

# Set default build type to "Release", change it in the GUI if you need to build with debug.
IF(NOT CMAKE_BUILD_TYPE)
	SET(CMAKE_BUILD_TYPE Release CACHE STRING
		"Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
	FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

#Some screen output
MESSAGE(STATUS "OS detected: ${CMAKE_SYSTEM_NAME}")
MESSAGE(STATUS "CXX Compiler detected: ${CMAKE_CXX_COMPILER_ID}")

# Setting OS specific cmake variables
# 1 - OSX
IF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	SET(CMAKE_MACOSX_RPATH OFF)
	SET(CMAKE_OSX_ROOT "")
	SET(CMAKE_OSX_DEPLOYMENT_TARGET "")
	# Homebrew directory added to path
	SET(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} /usr/local/lib)
ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
# 2 - NONE

# Some screen output
MESSAGE(STATUS "CMake additional search path for libraries: ${CMAKE_LIBRARY_PATH}")

# We set compiler specific flags
INCLUDE(cmake_scripts/set_compiler_flags.cmake)

# Provides build options to CMake.
# Build Option: when active headers will be installed in the CMAKE_INSTALL_PREFIX/include directory 
OPTION(INSTALL_HEADERS "Installs the header files" OFF)

# Build Option: when active the file main.cpp is built and linked to the PaGMO static library.
OPTION(BUILD_MAIN "Build 'main.cpp'." ON)

# Build Option: when active the bindings to Python are compiled and linked (installation).
OPTION(BUILD_PYGMO "Build Python bindings." OFF)

# Build Option: when active the GTOP database problems are built.
OPTION(ENABLE_GTOP_DATABASE "Build GTOP database problems (interplanetary transfers)." OFF)
IF(ENABLE_GTOP_DATABASE)
		ADD_DEFINITIONS(-DPAGMO_ENABLE_KEP_TOOLBOX)
ENDIF(ENABLE_GTOP_DATABASE)

# Build Option: minimisers from the GNU scientific library (GSL).
OPTION(ENABLE_GSL "Enable support for GSL minimisers (requires GSL >= 1.13)." OFF)

# Build Option: algorithms from the NLopt library.
OPTION(ENABLE_NLOPT "Enable support for NLopt minimisers." OFF)

# Build Option: SNOPT optimiser.
OPTION(ENABLE_SNOPT "Enable support for SNOPT minimiser." OFF)

# Build Option: IPOPT optimiser.
OPTION(ENABLE_IPOPT "Enable support for IPOPT minimiser." OFF)

# Build Option: WORHP optimiser.
OPTION(ENABLE_WORHP "Enable support for WORHP minimiser." OFF)

# Build Option: support for MPI clustering.
OPTION(ENABLE_MPI "Enable support for the Message Passage Interface (MPI)." OFF)

# Build option: enable test set.
OPTION(ENABLE_TESTS "Build test set." OFF)

# Build Option: build executable for the examples
OPTION(BUILD_EXAMPLES "Build examples." OFF)

SET(DYNAMIC_LIB_PAGMO_USE_FLAGS "-DBOOST_THREAD_USE_DLL -DBOOST_SERIALIZATION_DYN_LINK=1")
# NOTE: for system Boost, we are always going to use the system DLLs.
SET(STATIC_LIB_PAGMO_USE_FLAGS "-DBOOST_THREAD_USE_DLL -DBOOST_SERIALIZATION_DYN_LINK=1")
SET(STATIC_LIB_PAGMO_BUILD_FLAGS "-DBOOST_SERIALIZATION_DYN_LINK=1")
SET(DYNAMIC_LIB_PAGMO_BUILD_FLAGS "-DBOOST_SERIALIZATION_DYN_LINK=1")

# Setting the boost libraries needed for PaGMO
SET(REQUIRED_BOOST_LIBS system serialization thread)
IF(BUILD_PYGMO)
	INCLUDE(cmake_scripts/python_setup.cmake)
	IF(PYTHON_VERSION_MAJOR LESS 3)
		SET(REQUIRED_BOOST_LIBS ${REQUIRED_BOOST_LIBS} python)
	ELSE(PYTHON_VERSION_MAJOR LESS 3)
		SET(REQUIRED_BOOST_LIBS ${REQUIRED_BOOST_LIBS} python3)
	ENDIF(PYTHON_VERSION_MAJOR LESS 3)
ENDIF(BUILD_PYGMO)

IF(ENABLE_GTOP_DATABASE)
	SET(REQUIRED_BOOST_LIBS ${REQUIRED_BOOST_LIBS} date_time)
ENDIF(ENABLE_GTOP_DATABASE)
MESSAGE(STATUS "Required Boost libraries: ${REQUIRED_BOOST_LIBS}")

FIND_PACKAGE(Boost 1.48.0 REQUIRED COMPONENTS "${REQUIRED_BOOST_LIBS}")
MESSAGE(STATUS "Detected Boost version: ${Boost_VERSION}")
# Include system Boost headers.
MESSAGE(STATUS "Boost include dirs: ${Boost_INCLUDE_DIRS}")
MESSAGE(STATUS "Boost libraries: ${Boost_LIBRARIES}")
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
SET(MANDATORY_BOOST_LIBS ${Boost_SYSTEM_LIBRARY} ${Boost_THREAD_LIBRARY} ${Boost_SERIALIZATION_LIBRARY})
IF(ENABLE_GTOP_DATABASE)
	SET(MANDATORY_BOOST_LIBS ${MANDATORY_BOOST_LIBS} ${Boost_DATE_TIME_LIBRARY})
ENDIF(ENABLE_GTOP_DATABASE)
SET(MANDATORY_LIBRARIES ${MANDATORY_LIBRARIES} ${MANDATORY_BOOST_LIBS})

IF(ENABLE_GTOP_DATABASE)
	# We look for keplerian_toolbox library and headers
	FIND_LIBRARY(KEP_TOOLBOX_LIB NAMES keplerian_toolbox)
	FIND_PATH(KEP_TOOLBOX_HEADERS NAMES keplerian_toolbox/keplerian_toolbox.h)
	IF(NOT KEP_TOOLBOX_LIB OR NOT KEP_TOOLBOX_HEADERS)
		MESSAGE(FATAL_ERROR "I cannot find the keplerian_toolbox library or headers, please install it and make sure it can be found")
	ENDIF(NOT KEP_TOOLBOX_LIB OR NOT KEP_TOOLBOX_HEADERS)
	MESSAGE(STATUS "Keplerian Toolbox library: ${KEP_TOOLBOX_LIB}")
	MESSAGE(STATUS "Keplerian Toolbox include dir: ${KEP_TOOLBOX_HEADERS}")
	SET(MANDATORY_LIBRARIES ${MANDATORY_LIBRARIES} ${KEP_TOOLBOX_LIB})
	INCLUDE_DIRECTORIES(${KEP_TOOLBOX_HEADERS})
ENDIF(ENABLE_GTOP_DATABASE)

# Initial thread setup.
FIND_PACKAGE(Threads REQUIRED)
MESSAGE(STATUS "Thread library: ${CMAKE_THREAD_LIBS_INIT}")
SET(MANDATORY_LIBRARIES ${MANDATORY_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})

# Some common platform switches.
SET(PAGMO_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX})
IF(UNIX)
	# Install path for libraries.
	SET(LIB_INSTALL_PATH "lib")
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
	IF(MINGW)
		# Linker and compiler flags needed for the support of multithread applications.
		SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS_INIT} -mthreads")
		SET(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS_INIT} -mthreads")
		SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS_INIT} -mthreads")
		SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mthreads")
	ENDIF(MINGW)
	IF(MSVC)
		#This flags are necessary for MSVC requires boost libs to be named in a stupid way. (i.e. without lib in front)
		ADD_DEFINITIONS(-DBOOST_ALL_DYN_LINK)
		#This flag is necesary for MSVC to access mathematical constants such as M_PI,...
		ADD_DEFINITIONS(-D_USE_MATH_DEFINES)
		#This is necessary to add Additional Library Directories in the linker path of MSVC
		link_directories(${Boost_LIBRARY_DIRS})
	ENDIF(MSVC)
ENDIF(UNIX)



# BLAS setup for PaGMO.
SET(BLAS_IS_SETUP FALSE)
MACRO(SETUP_BLAS_FOR_PAGMO)
	IF(NOT BLAS_IS_SETUP)
		MESSAGE(STATUS "BLAS support was requested, identifying system BLAS libraries.")
		ENABLE_LANGUAGE(Fortran)
		FIND_PACKAGE(BLAS REQUIRED)
		MESSAGE(STATUS "BLAS libraries: ${BLAS_LIBRARIES}")
		MESSAGE(STATUS "BLAS linker flags: ${BLAS_LINKER_FLAGS}")
		SET(MANDATORY_LIBRARIES ${MANDATORY_LIBRARIES} ${BLAS_LIBRARIES})
		SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS_INIT} ${BLAS_LINKER_FLAGS}")
		SET(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS_INIT} ${BLAS_LINKER_FLAGS}")
		SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS_INIT} ${BLAS_LINKER_FLAGS}")
		SET(BLAS_IS_SETUP TRUE)
	ENDIF(NOT BLAS_IS_SETUP)
ENDMACRO(SETUP_BLAS_FOR_PAGMO)

IF(ENABLE_MPI)
	INCLUDE(FindMPI)
	IF(NOT MPI_FOUND)
	MESSAGE(FATAL_ERROR "MPI support was requested, but the MPI environment was not detected. Please check your setup.")
	ENDIF(NOT MPI_FOUND)
	MESSAGE(STATUS "MPI libraries: ${MPI_LIBRARIES}")
	MESSAGE(STATUS "MPI include path: ${MPI_INCLUDE_PATH}")
	MESSAGE(STATUS "MPI compile flags: ${MPI_COMPILE_FLAGS}")
	MESSAGE(STATUS "MPI link flags: ${MPI_LINK_FLAGS}")
	MESSAGE(STATUS "MPI exec file: ${MPIEXEC}")
	MESSAGE(STATUS "MPI exec num_processors flag: ${MPIEXEC_NUMPROC_FLAG}")
	MESSAGE(STATUS "MPI exec pre flags: ${MPIEXEC_PREFLAGS}")
	MESSAGE(STATUS "MPI exec post flags: ${MPIEXEC_POSTFLAGS}")
	SET(MANDATORY_LIBRARIES ${MANDATORY_LIBRARIES} ${MPI_LIBRARIES})
	INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
	SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS_INIT} ${MPI_LINK_FLAGS}")
	SET(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS_INIT} ${MPI_LINK_FLAGS}")
	SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS_INIT} ${MPI_LINK_FLAGS}")
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_COMPILE_FLAGS}")
	ADD_DEFINITIONS(-DPAGMO_ENABLE_MPI)
ENDIF(ENABLE_MPI)

# If GSL support is requested, look for the library.
IF(ENABLE_GSL)
	FIND_LIBRARY(GSL_GSL_LIBRARY NAMES gsl)
	FIND_LIBRARY(GSL_GSLCBLAS_LIBRARY NAMES gslcblas)
	FIND_PATH(GSL_INCLUDE_DIR NAMES gsl/gsl_version.h)
	MESSAGE(STATUS "GSL library: ${GSL_GSL_LIBRARY}")
	MESSAGE(STATUS "GSL CBLAS library: ${GSL_GSLCBLAS_LIBRARY}")
	MESSAGE(STATUS "GSL include dir: ${GSL_INCLUDE_DIR}")
	IF(NOT GSL_GSL_LIBRARY OR NOT GSL_INCLUDE_DIR)
		MESSAGE(FATAL_ERROR "GSL support was requested, but either the library or the include files could not be located.")
	ENDIF(NOT GSL_GSL_LIBRARY OR NOT GSL_INCLUDE_DIR)
	MESSAGE(STATUS "Running the GSL version checker.")
	TRY_RUN(SUFFICIENT_GSL_VERSION_RUN SUFFICIENT_GSL_VERSION_COMPILE ${CMAKE_BINARY_DIR}/compile_tests/ ${CMAKE_CURRENT_SOURCE_DIR}/cmake_scripts/gsl_version.cpp
		COMPILE_DEFINITIONS "${CMAKE_CCXX_FLAGS} -I\"${Boost_INCLUDE_DIRS}\" -I\"${GSL_INCLUDE_DIR}\"" RUN_OUTPUT_VARIABLE SUFFICIENT_GSL_VERSION_OUTPUT)
	MESSAGE(STATUS "GSL version checker output: ${SUFFICIENT_GSL_VERSION_OUTPUT}")
	IF(NOT SUFFICIENT_GSL_VERSION_COMPILE)
		MESSAGE(FATAL_ERROR "Error compiling the GSL version checker.")
	ENDIF(NOT SUFFICIENT_GSL_VERSION_COMPILE)
	IF(NOT ${SUFFICIENT_GSL_VERSION_RUN} EQUAL 0)
		MESSAGE(FATAL_ERROR "Error running the GSL version checker: either the GSL version is < 1.15 or your GSL installation is broken.")
	ELSE(NOT ${SUFFICIENT_GSL_VERSION_RUN} EQUAL 0)
		MESSAGE(STATUS "GSL version is fine.")
	ENDIF(NOT ${SUFFICIENT_GSL_VERSION_RUN} EQUAL 0)
	SET(MANDATORY_LIBRARIES ${MANDATORY_LIBRARIES} ${GSL_GSL_LIBRARY})
	IF(GSL_GSLCBLAS_LIBRARY)
		MESSAGE(STATUS "Found GSL BLAS library, linking it in.")
		SET(MANDATORY_LIBRARIES ${MANDATORY_LIBRARIES} ${GSL_GSLCBLAS_LIBRARY})
	ELSE(GSL_GSLCBLAS_LIBRARY)
		# If we could not locate GSL's own BLAS library, we will need an external one.
		MESSAGE(STATUS "Could not find GSL BLAS library, will need an external one.")
		SETUP_BLAS_FOR_PAGMO()
	ENDIF(GSL_GSLCBLAS_LIBRARY)
	ADD_DEFINITIONS(-DPAGMO_ENABLE_GSL)
	INCLUDE_DIRECTORIES("${GSL_INCLUDE_DIR}")
ENDIF(ENABLE_GSL)

# If NLopt support is requested, look for the library.
IF(ENABLE_NLOPT)
	FIND_LIBRARY(NLOPT_LIBRARY NAMES nlopt)
	FIND_PATH(NLOPT_INCLUDE_DIR NAMES nlopt.hpp)
	IF(NOT NLOPT_LIBRARY OR NOT NLOPT_INCLUDE_DIR)
		MESSAGE(FATAL_ERROR "NLopt support was requested, but the library could not be located. The library should be called 'libnlopt')")
	ENDIF(NOT NLOPT_LIBRARY OR NOT NLOPT_INCLUDE_DIR)
	MESSAGE(STATUS "NLopt library: ${NLOPT_LIBRARY}")
	MESSAGE(STATUS "NLopt include dir: ${NLOPT_INCLUDE_DIR}")
	SET(MANDATORY_LIBRARIES ${MANDATORY_LIBRARIES} ${NLOPT_LIBRARY})
	ADD_DEFINITIONS(-DPAGMO_ENABLE_NLOPT)
	INCLUDE_DIRECTORIES("${NLOPT_INCLUDE_DIR}")
ENDIF(ENABLE_NLOPT)

# If SNOPT support is requested, look for the libraries.
IF(ENABLE_SNOPT)
	SETUP_BLAS_FOR_PAGMO()

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

	FIND_LIBRARY(SNOPT_F2C_LIBRARY NAMES f2c)
	IF(NOT SNOPT_F2C_LIBRARY)
		MESSAGE(FATAL_ERROR "SNOPT support was requested, but the f2c library could not be located.")
	ENDIF(NOT SNOPT_F2C_LIBRARY)
	MESSAGE(STATUS "f2c library: ${SNOPT_F2C_LIBRARY}")

	FIND_LIBRARY(SNOPT_GFORTRAN_LIBRARY NAMES gfortran)
	IF(NOT SNOPT_GFORTRAN_LIBRARY)
		MESSAGE(FATAL_ERROR "SNOPT support was requested, but the gfortran library could not be located.")
	ENDIF(NOT SNOPT_GFORTRAN_LIBRARY)
	MESSAGE(STATUS "gfortran library: ${SNOPT_GFORTRAN_LIBRARY}")

		SET(MANDATORY_LIBRARIES ${MANDATORY_LIBRARIES} ${SNOPT_SNOPT_LIBRARY} ${SNOPT_SNPRINT_LIBRARY} ${SNOPT_F2C_LIBRARY} ${SNOPT_GFORTRAN_LIBRARY})

	ADD_DEFINITIONS(-DPAGMO_ENABLE_SNOPT)
ENDIF(ENABLE_SNOPT)

# If IPOPT support is requested, look for the libraries.
IF(ENABLE_IPOPT)
	#We look for ipopt libraries and headers
	FIND_LIBRARY(IPOPT_LIBRARY NAMES ipopt)
	FIND_PATH(IPOPT_INCLUDE_DIR NAMES coin/IpTNLP.hpp)
	MESSAGE(STATUS "IPOPT library: ${IPOPT_LIBRARY}")
	MESSAGE(STATUS "IPOPT include dir: ${IPOPT_INCLUDE_DIR}")
	IF(NOT IPOPT_LIBRARY OR NOT IPOPT_INCLUDE_DIR)
		MESSAGE(FATAL_ERROR "IPOPT support was requested, but ipopt library and headers could not be located. Please make sure that the linker can find the ipopt library and headers. The headers are assumed to be installed in a in coin/ directory")
	ENDIF(NOT IPOPT_LIBRARY OR NOT IPOPT_INCLUDE_DIR)

	IF(UNIX)
		#We set BLAS, LAPACK and -ldl
		SETUP_BLAS_FOR_PAGMO()

		#Note that LAPACK needs Fortran support, which has already been activated in SETUP_BLAS_FOR_PAGMO().
		FIND_PACKAGE(LAPACK REQUIRED)
		MESSAGE(STATUS "lapack libraries: ${LAPACK_LIBRARIES}")
		MESSAGE(STATUS "lapack linker flags: ${LAPACK_LINKER_FLAGS}")
		FIND_LIBRARY(DL_LIBRARY NAMES dl)
		MESSAGE(STATUS "dl library: ${DL_LIBRARY}")
		IF(NOT DL_LIBRARY)
			MESSAGE(FATAL_ERROR "IPOPT support was requested, but dl library could not be found.")
		ENDIF(NOT DL_LIBRARY)

		FIND_LIBRARY(COINHSL_LIBRARY NAMES coinhsl)

		FIND_LIBRARY(COINMETIS_LIBRARY NAMES coinmetis)

		IF(COINMETIS_LIBRARY)
			MESSAGE(STATUS "COINMETIS Library Found ... linking it in")
			MESSAGE(STATUS "coin metis library: ${COINMETIS_LIBRARY}")
			 SET(MANDATORY_LIBRARIES ${COINMETIS_LIBRARY} ${MANDATORY_LIBRARIES})
		ENDIF(COINMETIS_LIBRARY)

		FIND_LIBRARY(COINHSL_LIBRARY NAMES coinhsl)
		IF(COINHSL_LIBRARY)
			MESSAGE(STATUS "COINHSL Library Found ... linking it in")
			MESSAGE(STATUS "coin hsl library: ${COINHSL_LIBRARY}")
			SET(MANDATORY_LIBRARIES ${COINHSL_LIBRARY} ${MANDATORY_LIBRARIES})
		ENDIF(COINHSL_LIBRARY)

		FIND_LIBRARY(COINMUMPS_LIBRARY NAMES coinmumps)
		IF(COINMUMPS_LIBRARY)
			MESSAGE(STATUS "COINMUMPS Library Found ... linking it in")
			MESSAGE(STATUS "coin mumps library: ${COINMUMPS_LIBRARY}")
			SET(MANDATORY_LIBRARIES ${COINMUMPS_LIBRARY} ${MANDATORY_LIBRARIES})
		ENDIF(COINMUMPS_LIBRARY)

		IF(NOT COINHSL_LIBRARY AND NOT COINMUMPS_LIBRARY)
			MESSAGE(STATUS "Linear Solver for IPOPT could not be detected. I tried with HSL and MUMPS. You need to link manually the library if necessary")
		ENDIF(NOT COINHSL_LIBRARY AND NOT COINMUMPS_LIBRARY)
		SET(MANDATORY_LIBRARIES ${IPOPT_LIBRARY} ${MANDATORY_LIBRARIES} ${LAPACK_LIBRARIES} ${DL_LIBRARY})

	ELSE(UNIX)
		IF(NOT WIN32)
			MESSAGE(FATAL_ERROR "Sorry, at the moment pagmo can be built only on Unix and Win32 environments.")
		ENDIF(NOT WIN32)
		#If the platform is windows we just point to the ipopt precompiled dll wihch contains everything
		#The flag HAVE_CONFIG_H is activated as to include the config_ipopt.h header
		ADD_DEFINITIONS(-DHAVE_CONFIG_H)
		MESSAGE(STATUS "Windows is detected for ipopt linking you probably need to link manually to the libs")
		FIND_LIBRARY(COIN_BLAS_LIBRARY NAMES coinblas)
		FIND_LIBRARY(COIN_LAPACK_LIBRARY NAMES coinlapack)
		FIND_LIBRARY(COIN_MUMPS_LIBRARY NAMES coinmumps)
		SET(MANDATORY_LIBRARIES ${MANDATORY_LIBRARIES} ${IPOPT_LIBRARY} ${COIN_MUMPS_LIBRARY} ${COIN_BLAS_LIBRARY} ${COIN_LAPACK_LIBRARY} -lgfortran -lpthread)

	ENDIF(UNIX)

	ADD_DEFINITIONS(-DPAGMO_ENABLE_IPOPT)
	INCLUDE_DIRECTORIES("${IPOPT_INCLUDE_DIR}")
	SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS_INIT} ${LAPACK_LINKER_FLAGS}")
	SET(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS_INIT} ${LAPACK_LINKER_FLAGS}")
	SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS_INIT} ${LAPACK_LINKER_FLAGS}")
ENDIF(ENABLE_IPOPT)

# If WORHP support is requested, look for the libraries.
IF(ENABLE_WORHP)
	FIND_LIBRARY(WORHP_LIBRARY NAMES worhp)
	FIND_PATH(WORHP_INCLUDE_DIR NAMES worhp/worhp.h)
	IF(NOT WORHP_LIBRARY)
		MESSAGE(FATAL_ERROR "WORHP support was requested, but the library could not be located. The library should be called 'libworhp')")
	ENDIF(NOT WORHP_LIBRARY)
	IF(NOT WORHP_INCLUDE_DIR)
		MESSAGE(FATAL_ERROR "WORHP support was requested, but worhp/worhp.h not be located.")
	ENDIF(NOT WORHP_INCLUDE_DIR)
	MESSAGE(STATUS "WORHP library: ${WORHP_LIBRARY}")
	MESSAGE(STATUS "WORHP include dir: ${WORHP_INCLUDE_DIR}")
	INCLUDE_DIRECTORIES("${WORHP_INCLUDE_DIR}")
	SET(MANDATORY_LIBRARIES ${MANDATORY_LIBRARIES} ${WORHP_LIBRARY})
	ADD_DEFINITIONS(-DPAGMO_ENABLE_WORHP)
ENDIF(ENABLE_WORHP)


# Add the directory for the PaGMO library.
ADD_SUBDIRECTORY("${CMAKE_CURRENT_SOURCE_DIR}/src")

# Add the directory for the PyGMO library.
IF(BUILD_PYGMO)
	ADD_SUBDIRECTORY("${CMAKE_CURRENT_SOURCE_DIR}/PyGMO")
ENDIF(BUILD_PYGMO)

# From now on all targets will use the static PaGMO library. Add the relevant flags.
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${STATIC_LIB_PAGMO_USE_FLAGS}")

MESSAGE(STATUS "Build flags: " "${CMAKE_CXX_FLAGS}")
MESSAGE(STATUS "Module linker flags: " "${CMAKE_MODULE_LINKER_FLAGS}")
MESSAGE(STATUS "Shared linker flags: " "${CMAKE_SHARED_LINKER_FLAGS}")

# Link main to pagmo_static library.
IF(BUILD_MAIN)
	ADD_EXECUTABLE(main main.cpp)
		TARGET_LINK_LIBRARIES(main pagmo_static ${MANDATORY_LIBRARIES})
ENDIF(BUILD_MAIN)

IF(ENABLE_TESTS)
	ADD_SUBDIRECTORY("${CMAKE_CURRENT_SOURCE_DIR}/tests")
	FILE(COPY "${CMAKE_CURRENT_SOURCE_DIR}/tests/hypervolume_test_data/" DESTINATION "${CMAKE_BINARY_DIR}/tests/hypervolume_test_data/")
ENDIF(ENABLE_TESTS)

IF(BUILD_EXAMPLES)
	ADD_SUBDIRECTORY("${CMAKE_CURRENT_SOURCE_DIR}/examples")
ENDIF(BUILD_EXAMPLES)
