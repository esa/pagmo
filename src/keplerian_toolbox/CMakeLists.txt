PROJECT(Keplerian_Toolbox)

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.0)
# Setup and include CMake modules.
SET(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/../../cmake_modules)

# Configuration for GCC.
IF(CMAKE_COMPILER_IS_GNUCXX)
	# Add to the base flags extra warnings, and set up message reporting from the compiler so that it is Eclipse-friendly.
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wdisabled-optimization -fmessage-length=0")
	# Suggested for multithreaded code.
	ADD_DEFINITIONS(-D_REENTRANT)
ENDIF(CMAKE_COMPILER_IS_GNUCXX)

# Provides build options to CMake
# Build Option: when active the file main.cpp is built and linked to the PaGMO static library
OPTION(BUILD_MAIN "Build 'main.cpp'." ON)

# Let's include Boost's headers.
INCLUDE_DIRECTORIES("${CMAKE_SOURCE_DIR}/../../external_headers")

# Keplerian_toolbox lib source files.
SET(PAGMO_LIB_SRC_LIST
	${CMAKE_CURRENT_SOURCE_DIR}/epoch.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/planet.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/sims_flanagan/leg.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/sims_flanagan/fb_traj.cpp
	)

ADD_LIBRARY(keplerian_toolbox_static STATIC ${PAGMO_LIB_SRC_LIST})

# Link main to pagmo_static library.
IF(BUILD_MAIN)
	ADD_EXECUTABLE(main main.cpp)
	TARGET_LINK_LIBRARIES(main keplerian_toolbox_static)
ENDIF(BUILD_MAIN)

