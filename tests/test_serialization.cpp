/*****************************************************************************
 *   Copyright (C) 2004-2009 The PaGMO development team,                     *
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

#include <algorithm>
#include <vector>
#include "../src/pagmo.h"
#include "test_functions.h"
#include <string>
#include "../src/AstroToolbox/mga_dsm.h"
////////////
//new Libraries
#include <iostream>
#include <fstream>
// New Boost libararies
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include <numeric>

using namespace pagmo;


void save_de(const algorithm::de &alg, const char* filename)
{
    std::ofstream ofs( filename );
    boost::archive::text_oarchive oa( ofs );
    oa << alg;
}
void restore_de(algorithm::de & alg, const char* filename)
{
    std::ifstream ifs( filename );
    boost::archive::text_iarchive ia( ifs );
    ia >> alg;
}

void save_pso(const algorithm::pso &alg, const char* filename)
{
    std::ofstream ofs( filename );
    boost::archive::text_oarchive oa( ofs );
    oa << alg;
}
void restore_pso(algorithm::pso & alg, const char* filename)
{
    std::ifstream ifs( filename );
    boost::archive::text_iarchive ia( ifs );
    ia >> alg;
}

void save_cs(const algorithm::cs &alg, const char* filename)
{
    std::ofstream ofs( filename );
    boost::archive::text_oarchive oa( ofs );
    oa << alg;
}
void restore_cs(algorithm::cs & alg, const char* filename)
{
    std::ifstream ifs( filename );
    boost::archive::text_iarchive ia( ifs );
    ia >> alg;
}

void save_sga(const algorithm::sga &alg, const char* filename)
{
    std::ofstream ofs( filename );
    boost::archive::text_oarchive oa( ofs );
    oa << alg;
}
void restore_sga(algorithm::sga &alg, const char* filename)
{
    std::ifstream ifs( filename );
    boost::archive::text_iarchive ia( ifs );
    ia >> alg;
}
#ifdef PAGMO_ENABLE_IPOPT
void save_ipopt(const algorithm::ipopt &alg, const char* filename)
{
    std::ofstream ofs( filename );
    boost::archive::text_oarchive oa( ofs );
    oa << alg;
}
void restore_ipopt(algorithm::ipopt &alg, const char* filename)
{
    std::ifstream ifs( filename );
    boost::archive::text_iarchive ia( ifs );
    ia >> alg;
}
#endif

void save_schwefel(const problem::schwefel & prob, const char* filename)
{
    std::ofstream ofs( filename );
    boost::archive::text_oarchive oa( ofs );
    oa << prob;
}
void restore_schwefel(problem::schwefel & prob, const char* filename)
{
    std::ifstream ifs( filename );
    boost::archive::text_iarchive ia( ifs );
    ia >> prob;
}

void save_cassini_1(const problem::cassini_1 & prob, const char* filename)
{
    std::ofstream ofs( filename );
    boost::archive::text_oarchive oa( ofs );
    oa << prob;
}
void restore_cassini_1(problem::cassini_1 & prob, const char* filename)
{
    std::ifstream ifs( filename );
    boost::archive::text_iarchive ia( ifs );
    ia >> prob;
}

void save_luksan_vlcek_2(const problem::luksan_vlcek_2 & prob, const char* filename)
{
    std::ofstream ofs( filename );
    boost::archive::text_oarchive oa( ofs );
    oa << prob;
}
void restore_luksan_vlcek_2(problem::luksan_vlcek_2 & prob, const char* filename)
{
    std::ifstream ifs( filename );
    boost::archive::text_iarchive ia( ifs );
    ia >> prob;
}

void save_knapsack(const problem::knapsack & prob, const char* filename)
{
    std::ofstream ofs( filename );
    boost::archive::text_oarchive oa( ofs );
    oa << prob;
}
void restore_knapsack(problem::knapsack & prob, const char* filename)
{
    std::ifstream ifs( filename );
    boost::archive::text_iarchive ia( ifs );
    ia >> prob;
}

void save_tandem(const problem::tandem & prob, const char* filename)
{
    std::ofstream ofs( filename );
    boost::archive::text_oarchive oa( ofs );
    oa << prob;
}
void restore_tandem(problem::tandem & prob, const char* filename)
{
    std::ifstream ifs( filename );
    boost::archive::text_iarchive ia( ifs );
    ia >> prob;
}

void save_cassini_2(const problem::cassini_2 & prob, const char* filename)
{
    std::ofstream ofs( filename );
    boost::archive::text_oarchive oa( ofs );
    oa << prob;
}
void restore_cassini_2(problem::cassini_2 & prob, const char* filename)
{
    std::ifstream ifs( filename );
    boost::archive::text_iarchive ia( ifs );
    ia >> prob;
}

void save_population(const population & pop, const char* filename)
{
    std::ofstream ofs( filename );
    boost::archive::text_oarchive oa( ofs );
    oa << pop;
}
void restore_population(population & pop, const char* filename)
{
    std::ifstream ifs( filename );
    boost::archive::text_iarchive ia( ifs );
    ia >> pop;
}

void save_migration1(const migration::fair_r_policy & policy, const char* filename)
{
    std::ofstream ofs( filename );
    boost::archive::text_oarchive oa( ofs );
    oa << policy;
}
void restore_migration1(migration::fair_r_policy & policy, const char* filename)
{
    std::ifstream ifs( filename );
    boost::archive::text_iarchive ia( ifs );
    ia >> policy;
}

void save_migration2(const migration::worst_r_policy & policy, const char* filename)
{
    std::ofstream ofs( filename );
    boost::archive::text_oarchive oa( ofs );
    oa << policy;
}
void restore_migration2(migration::worst_r_policy & policy, const char* filename)
{
    std::ifstream ifs( filename );
    boost::archive::text_iarchive ia( ifs );
    ia >> policy;
}

void save_topology1(const topology::erdos_renyi & topo, const char* filename)
{
    std::ofstream ofs( filename );
    boost::archive::text_oarchive oa( ofs );
    oa << topo;
}
void restore_topology1(topology::erdos_renyi & topo, const char* filename)
{
    std::ifstream ifs( filename );
    boost::archive::text_iarchive ia( ifs );
    ia >> topo;
}

void save_topology2(const topology::watts_strogatz & topo, const char* filename)
{
    std::ofstream ofs( filename );
    boost::archive::text_oarchive oa( ofs );
    oa << topo;
}
void restore_topology2(topology::watts_strogatz & topo, const char* filename)
{
    std::ifstream ifs( filename );
    boost::archive::text_iarchive ia( ifs );
    ia >> topo;
}

int main()
{
	std::string fileName = "test.txt";
	/////
	// Testing some of the algorithms serialized
	/////

	// Algorithm DE
 	algorithm::de de_source = algorithm::de(500,.8,.9,2);
    algorithm::de de_dest = algorithm::de(600,.6,.7,1);
	std::cout << "Init. de algorithm: " << de_dest;
	//saving algorithm de1 into a text archive 
	save_de(de_source, fileName.c_str());
	//restoring the algoirhm de1 into a different object
	restore_de(de_dest, fileName.c_str());
	std::cout << "Rtrv. de algorithm: " << de_dest << endl << endl;

	// Algorithm PSO
    algorithm::pso pso_source = algorithm::pso(10,.65,2,2,.2,3);
    algorithm::pso pso_dest = algorithm::pso(15,.75,3,3,.4,2);
	std::cout << "Init. pso algorithm: " << pso_dest;
	//saving algorithm pso_source into a text archive 
	save_pso(pso_source, fileName.c_str());
	//restoring the algoirhm pso_source into a different object
	restore_pso(pso_dest, fileName.c_str());
	std::cout << "Rtrv. pso algorithm: " << pso_dest << endl << endl;

	// Algorithm CS
    algorithm::cs cs_source = algorithm::cs(1,.01,.1,.1);
    algorithm::cs cs_dest = algorithm::cs(99,.09);
	std::cout << "Init. cs algorithm: " << cs_dest;
	//saving algorithm cs_source into a text archive 
	save_cs(cs_source, fileName.c_str());
	//restoring the algoirhm cs_source into a different object
	restore_cs(cs_dest, fileName.c_str());
	std::cout << "Rtrv. cs algorithm: " << cs_dest << endl << endl;

	// Algorithm SGA
    algorithm::sga sga_source = algorithm::sga(10,.65,.2);
    algorithm::sga sga_dest = algorithm::sga(15,.75,.3);
	std::cout << "Init. sga algorithm: " << sga_dest;
	//saving algorithm sga_source into a text archive 
	save_sga(sga_source, fileName.c_str());
	//restoring the algoirhm sga_source into a different object
	restore_sga(sga_dest, fileName.c_str());
	std::cout << "Rtrv. sga algorithm: " << sga_dest << endl << endl;

#ifdef PAGMO_ENABLE_IPOPT
	// Algorithm IPOPT
	algorithm::ipopt ipopt_source = algorithm::ipopt(100);
    algorithm::ipopt ipopt_dest = algorithm::ipopt(5);
	std::cout << "Init. ipopt algorithm: " << ipopt_dest;
	//saving algorithm ipopt_source into a text archive 
	save_ipopt(ipopt_source_source, fileName.c_str());
	//restoring the algoirhm ipopt_source into a different object
	restore_ipopt(ipopt_dest, fileName.c_str());
	std::cout << "Rtrv. ipopt algorithm: " << ipopt_dest << endl << endl;
#endif


	/////
	// Testing some of the problems serialized
	/////

	problem::schwefel schwefel_source =  problem::schwefel(10);
	problem::schwefel schwefel_dest =  problem::schwefel(10);
	std::cout << "Inital schwefel problem: " << schwefel_dest << '\n';
	//saving algorithm de1 into a text archive 
	save_schwefel(schwefel_source, fileName.c_str());
	//restoring the algoirhm de1 into a different object
	restore_schwefel(schwefel_dest, fileName.c_str());
	std::cout << "schwefel problem retrieved form archive: " << schwefel_dest << '\n';
	
	problem::cassini_1 cassini_1_source =  problem::cassini_1();
	problem::cassini_1 cassini_1_dest =  problem::cassini_1();
	std::cout << "Inital cassini_1 problem: " << cassini_1_dest << '\n';
	//saving algorithm de1 into a text archive 
	save_cassini_1(cassini_1_source, fileName.c_str());
	//restoring the algoirhm de1 into a different object
	restore_cassini_1(cassini_1_dest, fileName.c_str());
	std::cout << "cassini_1 problem retrieved form archive: " << cassini_1_dest << '\n';

	problem::luksan_vlcek_2 luksan_vlcek_2_source =  problem::luksan_vlcek_2(16,1,5);
	problem::luksan_vlcek_2 luksan_vlcek_2_dest =  problem::luksan_vlcek_2(102,45,333);
	std::cout << "Inital luksan_vlcek_2 problem: " << luksan_vlcek_2_dest << '\n';
	//saving algorithm de1 into a text archive 
	save_luksan_vlcek_2(luksan_vlcek_2_source, fileName.c_str());
	//restoring the algoirhm de1 into a different object
	restore_luksan_vlcek_2(luksan_vlcek_2_dest, fileName.c_str());
	std::cout << "luksan_vlcek_2 problem retrieved form archive: " << luksan_vlcek_2_dest << '\n';


	double values_1[2], weights_1[2], values_2[3], weights_2[3];
	weights_1[0] = 1;
	weights_1[1] = 2;
	values_1[0] = 11;
	values_1[1] = 22;
	weights_2[2] = 2;
	weights_2[1] = 3;
	weights_2[0] = 4;
	values_2[2] = 22;
	values_2[1] = 33;
	values_2[0] = 65;
	problem::knapsack knapsack_source =  problem::knapsack(values_1,weights_1,6);
	problem::knapsack knapsack_dest =  problem::knapsack(values_2,weights_2,100);
	std::cout << "Initial knapsack problem: " << knapsack_dest << '\n';
	//saving algorithm de1 into a text archive 
	save_knapsack(knapsack_source, fileName.c_str());
	//restoring the algoirhm de1 into a different object
	restore_knapsack(knapsack_dest, fileName.c_str());
	std::cout << "knapsack problem retrieved form archive: " << knapsack_dest << '\n';

	problem::tandem tandem_source =  problem::tandem();
	problem::tandem tandem_dest =  problem::tandem();
	std::cout << "Initial tandem problem: " << tandem_dest << '\n';
	//saving algorithm de1 into a text archive 
	save_tandem(tandem_source, fileName.c_str());
	//restoring the algoirhm de1 into a different object
	restore_tandem(tandem_dest, fileName.c_str());
	std::cout << "tandem problem retrieved form archive: " << tandem_dest << '\n';

	problem::cassini_2 cassini_2_source =  problem::cassini_2();
	problem::cassini_2 cassini_2_dest =  problem::cassini_2();
	std::cout << "Inital cassini_2 problem: " << cassini_2_dest << '\n';
	//saving algorithm de1 into a text archive 
	save_cassini_2(cassini_2_source, fileName.c_str());
	//restoring the algoirhm de1 into a different object
	restore_cassini_2(cassini_2_dest, fileName.c_str());
	std::cout << "cassini_2 problem retrieved form archive: " << cassini_2_dest << '\n';

	// Testing the serialization of the population class
	population population_source = population(cassini_1_source,1);
	population population_dest = population(cassini_1_source,2);
	std::cout << "Inital population: " << population_dest << '\n';
	//saving algorithm de1 into a text archive 
	save_population(population_source, fileName.c_str());
	//restoring the algoirhm de1 into a different object
	restore_population(population_dest, fileName.c_str());
	std::cout << "cassini_2 problem retrieved form archive: " << population_dest << '\n';


	migration::fair_r_policy migration_source1 = migration::fair_r_policy();
	migration::fair_r_policy migration_dest1 = migration::fair_r_policy();
	std::cout << "Inital topology: " << migration_dest1 << '\n';
	//saving algorithm de1 into a text archive 
	save_migration1(migration_source1, fileName.c_str());
	//restoring the algoirhm de1 into a different object
	restore_migration1(migration_dest1, fileName.c_str());
	std::cout << "topology retrieved form archive: " << migration_dest1 << '\n';

	migration::worst_r_policy migration_source2 = migration::worst_r_policy();
	migration::worst_r_policy migration_dest2 = migration::worst_r_policy();
	std::cout << "Inital topology: " << migration_dest2 << '\n';
	//saving algorithm de1 into a text archive 
	save_migration2(migration_source2, fileName.c_str());
	//restoring the algoirhm de1 into a different object
	restore_migration2(migration_dest2, fileName.c_str());
	std::cout << "topology retrieved form archive: " << migration_dest2 << '\n';

	topology::erdos_renyi topology_source1 = topology::erdos_renyi(0.01);
	topology::erdos_renyi topology_dest1 = topology::erdos_renyi(0.02);
	std::cout << "Inital topology: " << topology_dest1 << '\n';
	//saving algorithm de1 into a text archive 
	save_topology1(topology_source1, fileName.c_str());
	//restoring the algoirhm de1 into a different object
	restore_topology1(topology_dest1, fileName.c_str());
	std::cout << "topology retrieved form archive: " << topology_dest1 << '\n';

	topology::watts_strogatz topology_source2 = topology::watts_strogatz(10,0.05,0);
	topology::watts_strogatz topology_dest2 = topology::watts_strogatz(10,0.03,0);
	std::cout << "Inital topology: " << topology_dest2 << '\n';
	//saving algorithm de1 into a text archive 
	save_topology2(topology_source2, fileName.c_str());
	//restoring the algoirhm de1 into a different object
	restore_topology2(topology_dest2, fileName.c_str());
	std::cout << "topology retrieved form archive: " << topology_dest2 << '\n';

	return 0;

}

