# Copyright (C) 2004-2009 The PaGMO development team,
# Advanced Concepts Team (ACT), European Space Agency (ESA)
# Copyright (c) 2014, K. Kumar (me@kartikkumar.com)
# Copyright (c) 2006, 2007 Montel Laurent, <montel@kde.org>
# Copyright (c) 2008, 2009 Gael Guennebaud, <g.gael@free.fr>
# Copyright (c) 2009 Benoit Jacob <jacob.benoit.1@gmail.com>
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
#
# This script tries to find the PaGMO library. This module supports requiring a minimum 
# version, e.g. you can do version, e.g. you can do find_package(PaGMO 1.1.5) to require
# version 1.1.5 or newer of PaGMO.
#
# Once done, this will define:
#
#     PAGMO_FOUND         - system has PaGMO lib with correct version;
#     PAGMO_INCLUDE_DIR   - the PaGMO include directory;
#     PAGMO_LIBRARY_DIR   - the PaGMO directory containing the library file(s);
#
# This script is based off of FindEigen3.cmake, which states that redistribution and use is 
# allowed according to the terms of the 2-clause BSD license.

macro(_pagmo_check_version)
  file(READ "${PAGMO_INCLUDE_DIR}/../VERSION" _pagmo_header)

  string(REGEX MATCH "define[ \t]+PAGMO_VERSION_MAJOR[ \t]+([0-9]+)" _pagmo_major_version_match "${_pagmo_header}")
  set(PAGMO_MAJOR_VERSION "${CMAKE_MATCH_1}")
  string(REGEX MATCH "define[ \t]+PAGMO_VERSION_MINOR[ \t]+([0-9]+)" _pagmo_minor_version_match "${_pagmo_header}")
  set(PAGMO_MINOR_VERSION "${CMAKE_MATCH_1}")  
  string(REGEX MATCH "define[ \t]+PAGMO_VERSION_PATCH[ \t]+([0-9]+)" _pagmo_patch_version_match "${_pagmo_header}")
  set(PAGMO_PATCH_VERSION "${CMAKE_MATCH_1}")  

  set(PAGMO_VERSION ${PAGMO_MAJOR_VERSION}.${PAGMO_MINOR_VERSION}.${PAGMO_PATCH_VERSION})
  if(${PAGMO_VERSION} VERSION_LESS ${PaGMO_FIND_VERSION})
    set(PAGMO_VERSION_OK FALSE)
  else(${PAGMO_VERSION} VERSION_LESS ${PaGMO_FIND_VERSION})
    set(PAGMO_VERSION_OK TRUE)
  endif(${PAGMO_VERSION} VERSION_LESS ${PaGMO_FIND_VERSION})

  if(NOT PAGMO_VERSION_OK)

    message(STATUS "PaGMO version ${PAGMO_VERSION} found in ${PAGMO_INCLUDE_DIR}, "
                   "but at least version ${PaGMO_FIND_VERSION} is required")
  endif(NOT PAGMO_VERSION_OK)

  set(PAGMO_LIBRARY "pagmo")
  link_directories(${PAGMO_LIBRARY_DIR})
endmacro(_pagmo_check_version)

if (PAGMO_INCLUDE_DIR)
 # in cache already
  _pagmo_check_version()
  set(PAGMO_FOUND ${PAGMO_VERSION_OK})

else (PAGMO_INCLUDE_DIR)
  find_path(PAGMO_BASE_PATH NAMES VERSION
      PATHS
      ${PROJECT_SOURCE_DIR}/../pagmo
      ${PROJECT_SOURCE_DIR}/../../pagmo
      ${PROJECT_SOURCE_DIR}/../../../pagmo
      ${PROJECT_SOURCE_DIR}/../../../../pagmo
      PATH_SUFFIXES pagmo 
    )
  set(PAGMO_INCLUDE_DIR ${PAGMO_BASE_PATH}/src/)
  set(PAGMO_LIBRARY_DIR ${PAGMO_BASE_PATH}/build/src/)

  if(PAGMO_INCLUDE_DIR)
    _pagmo_check_version()
  endif(PAGMO_INCLUDE_DIR)

  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(PaGMO DEFAULT_MSG PAGMO_INCLUDE_DIR PAGMO_VERSION_OK)

  mark_as_advanced(PAGMO_INCLUDE_DIR)

endif(PAGMO_INCLUDE_DIR)
