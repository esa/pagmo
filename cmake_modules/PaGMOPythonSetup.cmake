# Copyright (C) 2007, 2008 by Francesco Biscani
# bluescarni@gmail.com
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
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

INCLUDE(FindPythonLibs)

# Find Python libraries
FIND_PACKAGE(PythonLibs REQUIRED)
MESSAGE(STATUS "Python libraries: " "${PYTHON_LIBRARIES}")
INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})
MESSAGE(STATUS "Python library: " "${PYTHON_LIBRARY}")
SET(PYDEXTENSION FALSE)
SET(SOEXTENSION FALSE)
IF(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
	MESSAGE(STATUS "OSX system detected.")
	SET(PYTHON_MODULES_PATH "/Library/Python/2.5/site-packages")
	SET(SOEXTENSION TRUE)
ELSE(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
	IF(UNIX)
		STRING(REGEX MATCH libpython[0-9]\\.?[0-9] PYTHON_LIBRARY_VERSION_DOT ${PYTHON_LIBRARY})
		STRING(REGEX REPLACE libpython "" PYTHON_LIBRARY_VERSION_DOT ${PYTHON_LIBRARY_VERSION_DOT})
		STRING(REGEX REPLACE \\. "" PYTHON_LIBRARY_VERSION ${PYTHON_LIBRARY_VERSION_DOT})
		# In certain systems the Python lib is in /usr/lib(64), in others under /usr/lib/python2.5/config/.
		# We try here to catch both cases.
		#STRING(REGEX REPLACE "(python${PYTHON_LIBRARY_VERSION_DOT}/config/)?libpython.*"
		#	"python${PYTHON_LIBRARY_VERSION_DOT}/site-packages/" PYTHON_MODULES_PATH ${PYTHON_LIBRARY})
		# Let's use CMAKE_INSTALL_PREFIX, so that if we specify a different install path it will be respected.
		SET(PYTHON_MODULES_PATH ${CMAKE_INSTALL_PREFIX}/lib/python${PYTHON_LIBRARY_VERSION_DOT}/site-packages)
	ELSE(UNIX)
		STRING(REGEX MATCH python[0-9]\\.?[0-9] PYTHON_LIBRARY_VERSION ${PYTHON_LIBRARY})
		STRING(REGEX REPLACE python "" PYTHON_LIBRARY_VERSION ${PYTHON_LIBRARY_VERSION})
		SET(PYTHON_MODULES_PATH .)
	ENDIF(UNIX)
	MESSAGE(STATUS "Python library version: " ${PYTHON_LIBRARY_VERSION})
	IF(${PYTHON_LIBRARY_VERSION} GREATER 24 AND WIN32)
		MESSAGE(STATUS "Python >= 2.5 detected on WIN32 platform. Output extension for compiled modules will be '.pyd'.")
		SET(PYDEXTENSION TRUE)
	ENDIF(${PYTHON_LIBRARY_VERSION} GREATER 24 AND WIN32)
ENDIF(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
MESSAGE(STATUS "Python modules install path: " "${PYTHON_MODULES_PATH}")
