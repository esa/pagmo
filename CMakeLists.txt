PROJECT(PaGMO)

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.0)

# Add to the base flags extra warnings.
IF(CMAKE_COMPILER_IS_GNUCXX)
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -fmessage-length=0 -Wdisabled-optimization")
	# Suggested for multithreaded code.
	ADD_DEFINITIONS(-D_REENTRANT)
ENDIF(CMAKE_COMPILER_IS_GNUCXX)

# Setup and include CMake modules.
SET(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake_modules)
INCLUDE(CheckTypeSize)

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

OPTION(BUILD_MAIN "Build 'main.cpp'." ON)
OPTION(BUILD_MAIN_NG "Build 'main-ng.cpp'." ON)
OPTION(BUILD_PYGMO "Build Python bindings." OFF)

INCLUDE_DIRECTORIES("${CMAKE_SOURCE_DIR}")
INCLUDE_DIRECTORIES("${CMAKE_SOURCE_DIR}/src/AstroToolbox")
INCLUDE_DIRECTORIES("${CMAKE_SOURCE_DIR}/src/Functions/rng")
INCLUDE_DIRECTORIES("${CMAKE_SOURCE_DIR}/src/Functions/objfuns")
INCLUDE_DIRECTORIES("${CMAKE_SOURCE_DIR}/src/GOclasses/algorithms")
INCLUDE_DIRECTORIES("${CMAKE_SOURCE_DIR}/src/GOclasses/basic")
INCLUDE_DIRECTORIES("${CMAKE_SOURCE_DIR}/src/GOclasses/problems")
INCLUDE_DIRECTORIES("${CMAKE_SOURCE_DIR}/src/SolversThreads")

# Let's include Boost's headers.
INCLUDE_DIRECTORIES("${CMAKE_SOURCE_DIR}/external_headers")

SET(MANDATORY_LIBRARIES "")

# Some common platform switches.
SET(PAGMO_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX})
IF(UNIX)
	# Install path for libraries.
	SET(LIB_INSTALL_PATH "lib")
	# TODO: use FindThreads library from CMake here instead.
	# If we are in Unix, let's require pthread.
	FIND_LIBRARY(PTHREAD_LIBRARY pthread)
	IF(NOT PTHREAD_LIBRARY)
		MESSAGE(FATAL_ERROR "POSIX thread library not found, please install it or specify manually where it is located.")
	ENDIF(NOT PTHREAD_LIBRARY)
	SET(MANDATORY_LIBRARIES ${MANDATORY_LIBRARIES} pthread)
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
ELSE(UNIX)
	IF(NOT WIN32)
		MESSAGE(FATAL_ERROR "Sorry, at the moment pagmo can be built only on Unix and Win32 environments.")
	ENDIF(NOT WIN32)
	SET(LIB_INSTALL_PATH .)
	SET(PAGMO_INSTALL_PREFIX "${PAGMO_INSTALL_PREFIX} ${PAGMO_VERSION}")
	ADD_DEFINITIONS(-DPAGMO_WIN32)
	IF(MINGW)
		SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mthreads")
		SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -mthreads")
		SET(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} -mthreads")
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
ENDIF(BUILD_PYGMO)

# Link main to pagmo_static library.
IF(BUILD_MAIN)
	ADD_EXECUTABLE(main main.cpp)
	TARGET_LINK_LIBRARIES(main pagmo_static ${MANDATORY_LIBRARIES})
ENDIF(BUILD_MAIN)

# Link main-ng to pagmo_static library.
IF(BUILD_MAIN_NG)
	ADD_EXECUTABLE(main-ng main-ng.cpp)
	TARGET_LINK_LIBRARIES(main-ng pagmo_static ${MANDATORY_LIBRARIES})
ENDIF(BUILD_MAIN_NG)
