PROJECT(PaGMO)

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.0)
# Setup and include CMake modules.
SET(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake_modules)

INCLUDE(CheckTypeSize)
INCLUDE(PaGMOCPackSetup)

# Configuration for GCC.
IF(CMAKE_COMPILER_IS_GNUCXX)
	# Add to the base flags extra warnings, and set up message reporting from the compiler so that it is Eclipse-friendly.
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wdisabled-optimization -fmessage-length=0")
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

OPTION(BUILD_MAIN "Build 'main.cpp'." ON)
OPTION(BUILD_PYGMO "Build Python bindings." OFF)

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
	TARGET_LINK_LIBRARIES(main pagmo_static ${MANDATORY_LIBRARIES})
ENDIF(BUILD_MAIN)

INCLUDE(CPack)
