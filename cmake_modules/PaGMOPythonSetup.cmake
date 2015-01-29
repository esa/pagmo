# Copyright (C) 2004-2009 The PaGMO development team,
# Advanced Concepts Team (ACT), European Space Agency (ESA)
# http://apps.sourceforge.net/mediawiki/pagmo
# http://apps.sourceforge.net/mediawiki/pagmo/index.php?title=Developers
# http://apps.sourceforge.net/mediawiki/pagmo/index.php?title=Credits
# act@esa.int
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the
# Free Software Foundation, Inc.,
# 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.


# Only automate if finder mode is on.
IF(PYGMO_PYTHON_FINDER)
	# Only unset and refind if PYGMO_PYTHON_VARIABLE changed.
	IF(NOT ${PYGMO_PYTHON_VERSION_CACHED} STREQUAL ${PYGMO_PYTHON_VERSION})
		# Unset vars from cache
		UNSET(PYTHONINTERP_FOUND CACHE)
		UNSET(PYTHON_EXECUTABLE CACHE)
		UNSET(PYTHONLIBS_FOUND CACHE)
		UNSET(PYTHON_LIBRARIES CACHE)
		UNSET(PYTHON_INCLUDE_DIRS CACHE)
		UNSET(PYTHONLIBS_VERSION_STRING CACHE)

		# Deprecated variables
		UNSET(PYTHON_INCLUDE_PATH CACHE)
		UNSET(PYTHON_DEBUG_LIBRARIES CACHE)

		# Manual variables
		UNSET(PYTHON_LIBRARY CACHE)
		UNSET(PYTHON_INCLUDE_DIR CACHE)

		# Find Python executable and libraries
		FIND_PACKAGE(PythonInterp ${PYGMO_PYTHON_VERSION})
		FIND_PACKAGE(PythonLibs ${PYGMO_PYTHON_VERSION})

		# Check if version found matches version asked for
		IF(NOT ${PYTHON_VERSION_STRING} STREQUAL ${PYGMO_PYTHON_VERSION} OR NOT PYTHONINTERP_FOUND)
			UNSET(PYGMO_PYTHON_FINDER CACHE)
			SET(PYGMO_PYTHON_FINDER OFF CACHE BOOL "Find automatically Python interpreter and libraries based on PYGMO_PYTHON_VERSION")
			mark_as_advanced(PYGMO_PYTHON_VERSION_CACHED)
        		MESSAGE(FATAL_ERROR "Unable to locate requested Python ${PYGMO_PYTHON_VERSION}. Found ${PYTHON_VERSION_STRING} instead. \n\r 1) Choose a different version (${PYTHON_VERSION_STRING}) and re-enable PYGMO_PYTHON_FINDER.\n\r 2) Set PYTHON vars manually. Please make sure to set PYGMO_PYTHON_VERSION to match your Python major version (2 or 3) for Boost compatibility.")
		ENDIF(NOT ${PYTHON_VERSION_STRING} STREQUAL ${PYGMO_PYTHON_VERSION} OR NOT PYTHONINTERP_FOUND)
		
	ENDIF(NOT ${PYGMO_PYTHON_VERSION_CACHED} STREQUAL ${PYGMO_PYTHON_VERSION})

	# Set cached variable, so above is not rerun every time.
	UNSET(PYGMO_PYTHON_VERSION_CACHED CACHE)
	SET(PYGMO_PYTHON_VERSION_CACHED ${PYGMO_PYTHON_VERSION} CACHE STRING "Track changes to PYGMO_PYTHON_VERSION. DO NOT CHANGE.")
	mark_as_advanced(PYGMO_PYTHON_VERSION_CACHED)
ELSE(PYGMO_PYTHON_FINDER)
	# Bare minimum of variables for manual installation
	IF(NOT DEFINED PYTHON_LIBRARY)
		MESSAGE(FATAL_ERROR "Please set PYTHON_LIBRARY")
	ENDIF(NOT DEFINED PYTHON_LIBRARY)
	IF(NOT DEFINED PYTHON_INCLUDE_DIR)
		MESSAGE(FATAL_ERROR "Please set PYTHON_INCLUDE_DIR")
	ENDIF(NOT DEFINED PYTHON_INCLUDE_DIR)
ENDIF(PYGMO_PYTHON_FINDER)

# Need to run this every time
FIND_PACKAGE(PythonLibs)

# Give some feedback (values after / include deprecated CMake var
MESSAGE(STATUS "Python interpreter: " "${PYTHON_EXECUTABLE}")
MESSAGE(STATUS "Python library: " "${PYTHON_LIBRARIES}")
MESSAGE(STATUS "Python include path: " "${PYTHON_INCLUDE_DIRS}")

INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_DIRS})

# These flags are used to signal the need to override the default extension of the Python modules
# depending on the architecture. Under Windows, for instance, CMake produces shared objects as
# .dll files, but Python from 2.5 onwards requires .pyd files (hence the need to override).
# A similar thing happens in SuckOSX.
SET(PYDEXTENSION FALSE)
SET(SOEXTENSION FALSE)

IF(UNIX)

	# Sanity checks for manual pointing to python
	# (1) Python interpreter needed in order to detect the appropriate directory of modules installation.
	IF(NOT DEFINED PYTHON_EXECUTABLE)
		MESSAGE(FATAL_ERROR "Please set PYTHON_EXECUTABLE")
	ENDIF(NOT DEFINED PYTHON_EXECUTABLE)
	# (2) Python library should be set
	IF(NOT DEFINED PYTHON_LIBRARY)
		MESSAGE(FATAL_ERROR "Please set PYTHON_LIBRARY")
	ENDIF(NOT DEFINED PYTHON_LIBRARY)
	# (3) If libraries were never found PYTHON_LIBRARY_VERSION_DOT is unset
	IF(NOT DEFINED PYTHON_LIBRARY_VERSION_DOT)
		EXECUTE_PROCESS(COMMAND ${PYTHON_EXECUTABLE} -c "from distutils.sysconfig import get_python_lib; print(\"\".join([str(s) for s in get_python_lib().split(\"/\")[-2] if s.isdigit()]))" OUTPUT_VARIABLE PYTHON_LIBRARY_VERSION_DOT OUTPUT_STRIP_TRAILING_WHITESPACE)
	ENDIF(NOT DEFINED PYTHON_LIBRARY_VERSION_DOT)

	# Now we must establish if the installation dir for Python modules is named 'site-packages' (as usual)
	# or 'dist-packages' (apparently Ubuntu 9.04 or maybe Python 2.6, it's not clear).
	EXECUTE_PROCESS(COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_SOURCE_DIR}/cmake_modules/python_packages_dir.py
		OUTPUT_VARIABLE PY_PACKAGES_DIR OUTPUT_STRIP_TRAILING_WHITESPACE)
	MESSAGE(STATUS "Python packages dir is: ${PY_PACKAGES_DIR}")
	# SuckOSX suckages.
	IF(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
		MESSAGE(STATUS "OSX system detected.")
		SET(SOEXTENSION TRUE)
		# Let's determine Python version by running the interpreter with the --version flag.
                EXECUTE_PROCESS(COMMAND ${PYTHON_EXECUTABLE} --version OUTPUT_VARIABLE PY_VERSION_OSX ERROR_VARIABLE PY_VERSION_OSX OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_STRIP_TRAILING_WHITESPACE)
		MESSAGE(STATUS "Python interpeter returns string: " ${PY_VERSION_OSX})
		STRING(REGEX MATCH [0-9]*\\.[0-9]* PYTHON_LIBRARY_VERSION_DOT ${PY_VERSION_OSX})
	ELSE(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
		# In sane Unix system we can fetch the Python version number directly from the library.
		STRING(REGEX MATCH libpython[0-9]*\\.[0-9]* PYTHON_LIBRARY_VERSION_DOT ${PYTHON_LIBRARY})
	ENDIF(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
	# Remove the dot from the Python version.
	STRING(REGEX REPLACE libpython "" PYTHON_LIBRARY_VERSION_DOT ${PYTHON_LIBRARY_VERSION_DOT})
	STRING(REGEX REPLACE \\. "" PYTHON_LIBRARY_VERSION ${PYTHON_LIBRARY_VERSION_DOT})
	# Let's use CMAKE_INSTALL_PREFIX, so that if we specify a different install path it will be respected.
	SET(PYTHON_MODULES_PATH lib/python${PYTHON_LIBRARY_VERSION_DOT}/${PY_PACKAGES_DIR})
ELSE(UNIX)
	STRING(REGEX MATCH python[0-9]* PYTHON_LIBRARY_VERSION ${PYTHON_LIBRARY})
	STRING(REGEX REPLACE python "" PYTHON_LIBRARY_VERSION ${PYTHON_LIBRARY_VERSION})
	SET(PYTHON_MODULES_PATH .)
	IF(${PYTHON_LIBRARY_VERSION} GREATER 24 AND WIN32)
		MESSAGE(STATUS "Python >= 2.5 detected on WIN32 platform. Output extension for compiled modules will be '.pyd'.")
		SET(PYDEXTENSION TRUE)
	ENDIF(${PYTHON_LIBRARY_VERSION} GREATER 24 AND WIN32)
ENDIF(UNIX)
MESSAGE(STATUS "Python library version: " ${PYTHON_LIBRARY_VERSION})
MESSAGE(STATUS "Python modules install path: " "${PYTHON_MODULES_PATH}")
