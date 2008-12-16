PROJECT(PaGMO)

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.0)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -fmessage-length=0 -Wdisabled-optimization")
SET(CMAKE_BUILD_TYPE "Release")

INCLUDE(CheckTypeSize)

INCLUDE_DIRECTORIES("${CMAKE_SOURCE_DIR}")
INCLUDE_DIRECTORIES("${CMAKE_SOURCE_DIR}/AstroToolbox")
INCLUDE_DIRECTORIES("${CMAKE_SOURCE_DIR}/Functions/rng")
INCLUDE_DIRECTORIES("${CMAKE_SOURCE_DIR}/Functions/objfuns")
INCLUDE_DIRECTORIES("${CMAKE_SOURCE_DIR}/GOclasses/algorithms")
INCLUDE_DIRECTORIES("${CMAKE_SOURCE_DIR}/GOclasses/basic")
INCLUDE_DIRECTORIES("${CMAKE_SOURCE_DIR}/GOclasses/problems")
INCLUDE_DIRECTORIES("${CMAKE_SOURCE_DIR}/SolversThreads")

# Let's include Boost's headers.
INCLUDE_DIRECTORIES("${CMAKE_SOURCE_DIR}/external_headers")

SET(MANDATORY_LIBRARIES "")

CHECK_TYPE_SIZE("void *" POINTER_SIZE)
MESSAGE(STATUS "Pointer size = ${POINTER_SIZE}")

SET(PAGMO_LIB_SRC_LIST
	${CMAKE_SOURCE_DIR}/AstroToolbox/Astro_Functions.cpp
	${CMAKE_SOURCE_DIR}/AstroToolbox/Lambert.cpp
	${CMAKE_SOURCE_DIR}/AstroToolbox/mga.cpp
	${CMAKE_SOURCE_DIR}/AstroToolbox/mga_dsm.cpp
	${CMAKE_SOURCE_DIR}/AstroToolbox/misc4Tandem.cpp
	${CMAKE_SOURCE_DIR}/AstroToolbox/Pl_Eph_An.cpp
	${CMAKE_SOURCE_DIR}/AstroToolbox/PowSwingByInv.cpp
	${CMAKE_SOURCE_DIR}/AstroToolbox/propagateKEP.cpp
	${CMAKE_SOURCE_DIR}/AstroToolbox/time2distance.cpp
	${CMAKE_SOURCE_DIR}/Functions/rng/PkRandom.cpp
	${CMAKE_SOURCE_DIR}/Functions/objfuns/trajobjfuns.cpp
	${CMAKE_SOURCE_DIR}/GOclasses/algorithms/ASA.cpp
	${CMAKE_SOURCE_DIR}/GOclasses/algorithms/DE.cpp
	${CMAKE_SOURCE_DIR}/GOclasses/algorithms/MPSO.cpp
	${CMAKE_SOURCE_DIR}/GOclasses/algorithms/PSO.cpp
	${CMAKE_SOURCE_DIR}/GOclasses/algorithms/SGA.cpp
	${CMAKE_SOURCE_DIR}/GOclasses/basic/individual.cpp
	${CMAKE_SOURCE_DIR}/GOclasses/basic/population.cpp
	${CMAKE_SOURCE_DIR}/GOclasses/problems/GOproblem.cpp
	${CMAKE_SOURCE_DIR}/GOclasses/problems/TrajectoryProblems.cpp
	${CMAKE_SOURCE_DIR}/SolversThreads/SolversThreads.cpp
)

# Boost thread setup.
ADD_DEFINITIONS(-DBOOST_THREAD_BUILD_LIB)
SET(PAGMO_LIB_SRC_LIST ${PAGMO_LIB_SRC_LIST} boost_thread/tss_null.cpp)
IF(WIN32)
	SET(PAGMO_LIB_SRC_LIST ${PAGMO_LIB_SRC_LIST}
		boost_thread/win32/exceptions.cpp
		boost_thread/win32/thread.cpp
		boost_thread/win32/tss_dll.cpp
		boost_thread/win32/tss_pe.cpp
	)
ELSE(WIN32)
    # If we are not in Windows, let's require pthread.
    FIND_LIBRARY(PTHREAD_LIBRARY pthread)
    IF(NOT PTHREAD_LIBRARY)
        MESSAGE(FATAL_ERROR "pthread library not found, please specify manually where it is located.")
    ENDIF(NOT PTHREAD_LIBRARY)
    # TODO: check about the dl library here.
    SET(MANDATORY_LIBRARIES ${MANDATORY_LIBRARIES} pthread)
	SET(PAGMO_LIB_SRC_LIST ${PAGMO_LIB_SRC_LIST}
		boost_thread/pthread/exceptions.cpp
		boost_thread/pthread/once.cpp
		boost_thread/pthread/thread.cpp
	)
ENDIF(WIN32)

ADD_LIBRARY(pagmo STATIC ${PAGMO_LIB_SRC_LIST})

ADD_EXECUTABLE(main main.cpp)

TARGET_LINK_LIBRARIES(main pagmo ${MANDATORY_LIBRARIES})

