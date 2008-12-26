PROJECT(PaGMO)

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.0)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -fmessage-length=0 -Wdisabled-optimization")

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

# Suggested for multithreaded code.
ADD_DEFINITIONS(-D_REENTRANT)

SET(MANDATORY_LIBRARIES "")

# Threads setup.
IF(NOT WIN32)
	# If we are not in Windows, let's require pthread.
	FIND_LIBRARY(PTHREAD_LIBRARY pthread)
	IF(NOT PTHREAD_LIBRARY)
		MESSAGE(FATAL_ERROR "POSIX thread library not found, please install it or specify manually where it is located.")
	ENDIF(NOT PTHREAD_LIBRARY)
	# TODO: check about the dl library here.
	SET(MANDATORY_LIBRARIES ${MANDATORY_LIBRARIES} pthread)
ENDIF(NOT WIN32)

CHECK_TYPE_SIZE("void *" POINTER_SIZE)
MESSAGE(STATUS "Pointer size = ${POINTER_SIZE}")

ADD_SUBDIRECTORY("${CMAKE_SOURCE_DIR}/src")

# Link to pagmo_static library if main is requested.
IF(BUILD_MAIN)
	ADD_EXECUTABLE(main main.cpp)
	TARGET_LINK_LIBRARIES(main pagmo_static ${MANDATORY_LIBRARIES})
ENDIF(BUILD_MAIN)
