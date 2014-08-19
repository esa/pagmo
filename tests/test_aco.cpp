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
    // create an instance of aco
    pagmo::algorithm::aco simple(300, 150, 0.2);
    
    /** example of conversion from list to chromosome
     * 
     *  | 0 1 2
     * -|------
     * 0| - 0 1
     * 1| 1 - 0
     * 2| 0 1 -
     * 
     * 0 2 1 => -01 1-0 01-
     */        
    
    std::vector<size_t> list = {0, 2, 1};
    decision_vector dec_target = {0,1, 1,0, 0,1};
    
    if(dec_target != simple.tour2chromosome(list)) {
        std::cout << "the conversion from list to chromosome is incorrect!\n";
        return 1;
    }
    
    // create a tsp problem matrix (symmetric)
    std::vector<std::vector<double> > burma14 = {
        {    0.0,	153.0,	510.0,	 706.0,	  966.0,	581.0,	455.0,	 70.0,	160.0,	 372.0,	157.0,	567.0,	342.0,	398.0},
        {  153.0,	  0.0,	422.0,	 664.0,	  997.0,	598.0,	507.0,	197.0,	311.0,	 479.0,	310.0,	581.0,	417.0,	376.0},
        {  510.0,	422.0,	  0.0,	 289.0,	  744.0,	390.0,	437.0,	491.0,	645.0,	 880.0,	618.0,	374.0,	455.0,	211.0},
        {  706.0,	664.0,	289.0,	   0.0,	  491.0,	265.0,	410.0,	664.0,	804.0,  1070.0,	768.0,	259.0,	499.0,	310.0},
        {  966.0,	997.0,	744.0,	 491.0,	    0.0,	400.0,	514.0,	902.0,	990.0,  1261.0,	947.0,	418.0,	635.0,	636.0},
        {  581.0,	598.0,	390.0,	 265.0,	  400.0,	  0.0,	168.0,	522.0,	634.0,	 910.0,	593.0,	 19.0,	284.0,	239.0},
        {  455.0,	507.0,	437.0,	 410.0,	  514.0,	168.0,	  0.0,	389.0,	482.0,	 757.0,	439.0,	163.0,	124.0,	232.0},
        {   70.0,	197.0,	491.0,	 664.0,	  902.0,	522.0,	389.0,	  0.0,	154.0,	 406.0,	133.0,	508.0,	273.0,	355.0},
        {  160.0,	311.0,	645.0,	 804.0,	  990.0,	634.0,	482.0,	154.0,	  0.0,	 276.0,	 43.0,	623.0,	358.0,	498.0},
        {  372.0,	479.0,	880.0,  1070.0,  1261.0,        910.0,	757.0,	406.0,	276.0,	   0.0,	318.0,	898.0,	633.0,	761.0},
        {  157.0,	310.0,	618.0,	 768.0,	  947.0,	593.0,	439.0,	133.0,	 43.0,	 318.0,	  0.0,	582.0,	315.0,	464.0},
        {  567.0,	581.0,	374.0,	 259.0,	  418.0,	 19.0,	163.0,	508.0,	623.0,	 898.0,	582.0,	  0.0,	275.0,	221.0},
        {  342.0,	417.0,	455.0,	 499.0,	  635.0,	284.0,	124.0,	273.0,	358.0,	 633.0,	315.0,	275.0,	  0.0,	247.0},
        {  398.0,	376.0,	211.0,	 310.0,	  636.0,	239.0,	232.0,	355.0,	498.0,	 761.0,	464.0,	221.0,	247.0,	  0.0}};
    
    // create a tsp problem
    pagmo::problem::tsp prob(burma14);
    
    // create a population
    pagmo::population pop(prob, 100);
    
    // run the evolution method
    simple.evolve(pop);
    
    // print the champion
    std::cout << pop.champion();
    
    /*
from PyGMO.problem import tsp
from PyGMO.algorithm import aco
from PyGMO import population
from PyGMO.util import tsp as tsputil
xml = tsputil.read_tsplib('ulysses22.xml')
prob = tsp(xml)
pop = population(prob, 100)
algo = aco(90, 150)
result = algo.evolve(pop)
print result.champion
    */
    
    // all iz well
    return 0;
}