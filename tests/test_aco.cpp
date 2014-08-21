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
    pagmo::algorithm::aco simple(1500, 250, 0.5);
    
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

    // create a tsp problem matrix (symmetric) ulysses22
    std::vector<std::vector<double> > ulysses22 = {
{    0.0, 509.0, 501.0, 312.0,1019.0, 736.0, 656.0,  60.0,1039.0, 726.0,2314.0, 479.0, 448.0, 479.0,  619.0, 150.0, 342.0, 323.0, 635.0, 604.0, 596.0, 202.0 },
{  509.0,   0.0, 126.0, 474.0,1526.0,1226.0,1133.0, 532.0,1449.0,1122.0,2789.0, 958.0, 941.0, 978.0, 1127.0, 542.0, 246.0, 510.0,1047.0,1021.0,1010.0, 364.0 },
{  501.0, 126.0,   0.0, 541.0,1516.0,1184.0,1084.0, 536.0,1371.0,1045.0,2728.0, 913.0, 904.0, 946.0, 1115.0, 499.0, 321.0, 577.0, 976.0, 952.0, 941.0, 401.0 },
{  312.0, 474.0, 541.0,   0.0,1157.0, 980.0, 919.0, 271.0,1333.0,1029.0,2553.0, 751.0, 704.0, 720.0,  783.0, 455.0, 228.0,  37.0, 936.0, 904.0, 898.0, 171.0 },
{ 1019.0,1526.0,1516.0,1157.0,   0.0, 478.0, 583.0, 996.0, 858.0, 855.0,1504.0, 677.0, 651.0, 600.0,  401.0,1033.0,1325.0,1134.0, 818.0, 808.0, 820.0,1179.0 },
{  736.0,1226.0,1184.0, 980.0, 478.0,   0.0, 115.0, 740.0, 470.0, 379.0,1581.0, 271.0, 289.0, 261.0,  308.0, 687.0,1077.0, 970.0, 342.0, 336.0, 348.0, 932.0 },
{  656.0,1133.0,1084.0, 919.0, 583.0, 115.0,   0.0, 667.0, 455.0, 288.0,1661.0, 177.0, 216.0, 207.0,  343.0, 592.0, 997.0, 913.0, 236.0, 226.0, 237.0, 856.0 },
{   60.0, 532.0, 536.0, 271.0, 996.0, 740.0, 667.0,   0.0,1066.0, 759.0,2320.0, 493.0, 454.0, 479.0,  598.0, 206.0, 341.0, 278.0, 666.0, 634.0, 628.0, 194.0 },
{ 1039.0,1449.0,1371.0,1333.0, 858.0, 470.0, 455.0,1066.0,   0.0, 328.0,1387.0, 591.0, 650.0, 656.0,  776.0, 933.0,1367.0,1333.0, 408.0, 438.0, 447.0,1239.0 },
{  726.0,1122.0,1045.0,1029.0, 855.0, 379.0, 288.0, 759.0, 328.0,   0.0,1697.0, 333.0, 400.0, 427.0,  622.0, 610.0,1046.0,1033.0,  96.0, 128.0, 133.0, 922.0 },
{ 2314.0,2789.0,2728.0,2553.0,1504.0,1581.0,1661.0,2320.0,1387.0,1697.0,   0.0,1838.0,1868.0,1841.0, 1789.0,2248.0,2656.0,2540.0,1755.0,1777.0,1789.0,2512.0 },
{  479.0, 958.0, 913.0, 751.0, 677.0, 271.0, 177.0, 493.0, 591.0, 333.0,1838.0,   0.0,  68.0, 105.0,  336.0, 417.0, 821.0, 748.0, 243.0, 214.0, 217.0, 680.0 },
{  448.0, 941.0, 904.0, 704.0, 651.0, 289.0, 216.0, 454.0, 650.0, 400.0,1868.0,  68.0,   0.0,  52.0,  287.0, 406.0, 789.0, 698.0, 311.0, 281.0, 283.0, 645.0 },
{  479.0, 978.0, 946.0, 720.0, 600.0, 261.0, 207.0, 479.0, 656.0, 427.0,1841.0, 105.0,  52.0,   0.0,  237.0, 449.0, 818.0, 712.0, 341.0, 314.0, 318.0, 672.0 },
{  619.0,1127.0,1115.0, 783.0, 401.0, 308.0, 343.0, 598.0, 776.0, 622.0,1789.0, 336.0, 287.0, 237.0,    0.0, 636.0, 932.0, 764.0, 550.0, 528.0, 535.0, 785.0 },
{  150.0, 542.0, 499.0, 455.0,1033.0, 687.0, 592.0, 206.0, 933.0, 610.0,2248.0, 417.0, 406.0, 449.0,  636.0,   0.0, 436.0, 470.0, 525.0, 496.0, 486.0, 319.0 },
{  342.0, 246.0, 321.0, 228.0,1325.0,1077.0, 997.0, 341.0,1367.0,1046.0,2656.0, 821.0, 789.0, 818.0,  932.0, 436.0,   0.0, 265.0, 959.0, 930.0, 921.0, 148.0 },
{  323.0, 510.0, 577.0,  37.0,1134.0, 970.0, 913.0, 278.0,1333.0,1033.0,2540.0, 748.0, 698.0, 712.0,  764.0, 470.0, 265.0,   0.0, 939.0, 907.0, 901.0, 201.0 },
{  635.0,1047.0, 976.0, 936.0, 818.0, 342.0, 236.0, 666.0, 408.0,  96.0,1755.0, 243.0, 311.0, 341.0,  550.0, 525.0, 959.0, 939.0,   0.0,  33.0,  39.0, 833.0 },
{  604.0,1021.0, 952.0, 904.0, 808.0, 336.0, 226.0, 634.0, 438.0, 128.0,1777.0, 214.0, 281.0, 314.0,  528.0, 496.0, 930.0, 907.0,  33.0,   0.0,  14.0, 803.0 },
{  596.0,1010.0, 941.0, 898.0, 820.0, 348.0, 237.0, 628.0, 447.0, 133.0,1789.0, 217.0, 283.0, 318.0,  535.0, 486.0, 921.0, 901.0,  39.0,  14.0,   0.0, 794.0 },
{  202.0, 364.0, 401.0, 171.0,1179.0, 932.0, 856.0, 194.0,1239.0, 922.0,2512.0, 680.0, 645.0, 672.0,  785.0, 319.0, 148.0, 201.0, 833.0, 803.0, 794.0,   0.0 }};
    
    // create a tsp problem
    pagmo::problem::tsp prob(ulysses22);
    
    // create a population
    pagmo::population pop(prob, 200);
    
    // run the evolution method
    simple.evolve(pop);
    
    // print the champion
    std::cout << pop.champion().f << std::endl;
//    if(pop.champion().f[0] != 7013)
//        return 1;
    
    /*
     * 
     * NAME : ulysses22.opt.tour 
     * TYPE : TOUR
     * COMMENT : Optimal solution of ulysses22 (7013)
     * DIMENSION : 22
     * TOUR_SECTION 
     * 
     * 1, 14, 13, 12, 7, 6, 15, 5, 11, 9, 10, 19, 20, 21, 16, 3, 2, 17, 22, 4, 18, 8 
      
     * [0, 15, 12, 13, 11, 6, 5, 14, 4, 19, 20, 18, 9, 8, 10, 2, 1, 16, 21, 3, 17, 7]
     
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