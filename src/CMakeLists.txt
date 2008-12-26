SET(PAGMO_LIB_SRC_LIST
	${CMAKE_CURRENT_SOURCE_DIR}/AstroToolbox/Astro_Functions.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/AstroToolbox/Lambert.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/AstroToolbox/mga.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/AstroToolbox/mga_dsm.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/AstroToolbox/misc4Tandem.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/AstroToolbox/Pl_Eph_An.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/AstroToolbox/PowSwingByInv.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/AstroToolbox/propagateKEP.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/AstroToolbox/time2distance.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/Functions/objfuns/trajobjfuns.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/GOclasses/algorithms/ASA.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/GOclasses/algorithms/DE.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/GOclasses/algorithms/MPSO.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/GOclasses/algorithms/PSO.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/GOclasses/algorithms/SGA.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/GOclasses/basic/individual.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/GOclasses/basic/population.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/GOclasses/problems/TrajectoryProblems.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/SolversThreads/SolversThreads.cpp
)

# Boost thread setup.
ADD_DEFINITIONS(-DBOOST_THREAD_BUILD_LIB)
SET(PAGMO_LIB_SRC_LIST ${PAGMO_LIB_SRC_LIST} boost_thread/tss_null.cpp)
IF(WIN32)
	SET(PAGMO_LIB_SRC_LIST ${PAGMO_LIB_SRC_LIST}
		${CMAKE_CURRENT_SOURCE_DIR}/boost_thread/win32/exceptions.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/boost_thread/win32/thread.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/boost_thread/win32/tss_dll.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/boost_thread/win32/tss_pe.cpp
	)
ELSE(WIN32)
	SET(PAGMO_LIB_SRC_LIST ${PAGMO_LIB_SRC_LIST}
		${CMAKE_CURRENT_SOURCE_DIR}/boost_thread/pthread/exceptions.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/boost_thread/pthread/once.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/boost_thread/pthread/thread.cpp
	)
ENDIF(WIN32)

# Create a pagmo_static library if main is requested.
IF(BUILD_MAIN)
	ADD_LIBRARY(pagmo_static STATIC ${PAGMO_LIB_SRC_LIST})
ENDIF(BUILD_MAIN)
