PROJECT(PaGMO)
 
CMAKE_MINIMUM_REQUIRED(VERSION 2.4.0)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -fmessage-length=0 -Wdisabled-optimization")
SET(CMAKE_BUILD_TYPE "Release")

INCLUDE_DIRECTORIES("${CMAKE_SOURCE_DIR}")
INCLUDE_DIRECTORIES("${CMAKE_SOURCE_DIR}/AstroToolbox")
INCLUDE_DIRECTORIES("${CMAKE_SOURCE_DIR}/Functions/rng")
INCLUDE_DIRECTORIES("${CMAKE_SOURCE_DIR}/Functions/objfuns")
INCLUDE_DIRECTORIES("${CMAKE_SOURCE_DIR}/GOclasses/algorithms")
INCLUDE_DIRECTORIES("${CMAKE_SOURCE_DIR}/GOclasses/basic")
INCLUDE_DIRECTORIES("${CMAKE_SOURCE_DIR}/GOclasses/problems")
INCLUDE_DIRECTORIES("${CMAKE_SOURCE_DIR}/SolversThreads")

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
	${CMAKE_SOURCE_DIR}/AstroToolbox/ZeroFinder.cpp
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

ADD_LIBRARY(pagmo STATIC ${PAGMO_LIB_SRC_LIST})

ADD_EXECUTABLE(main main.cpp)

TARGET_LINK_LIBRARIES(main pagmo pthread)

