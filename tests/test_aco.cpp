 /*****************************************************************************
 *   Copyright (C) 2004-2013 The PaGMO development team,                     *
 *   Advanced Concepts Team (ACT), European Space Agency (ESA)               *
 *   http://apps.sourceforge.net/mediawiki/pagmo                             *
 *   http://apps.sourceforge.net/mediawiki/pagmo/index.php?title=Developers  *
 *   http://apps.sourceforge.net/mediawiki/pagmo/index.php?title=Credits     *
 *   act@esa.int                                                             *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program; if not, write to the                           *
 *   Free Software Foundation, Inc.,                                         *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.               *
 *****************************************************************************/
#include <iostream>
#include <iomanip>

#include "../src/pagmo.h"

using namespace pagmo;

int main()
{
    //if (test_simple() return 1;
    
    pagmo::algorithm::aco simple(100, 100, 0.2);
    
    /**
     *  | 0 1 2 3
     * -|--------
     * 0| 0 0 1 0
     * 1| 0 0 0 1
     * 2| 0 1 0 0
     * 3| 1 0 0 0
     * 
     */        
    
    std::vector<size_t> list = {0, 2, 1, 3};
    
    decision_vector dec_orig = {};
    decision_vector target = simple.list2decision_vector(list);
    
    //if(dec_orig != target) return 1;
    
    // all iz well
    return 0;
}