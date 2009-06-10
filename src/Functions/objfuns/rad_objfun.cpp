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

#include <cmath>

#include "rad_objfun.h"

using namespace std;

/*Continuous radiator function from Jose study*/

extern "C"
  {
    double real_radiator_extern_(double *, const int &);
  }

double radiator( const vector<double> &x )
{
        const int n = x.size();
        double J_result=0;
        const int number_layers=C_NUMBER_LAYERS;
        double* array = new double[number_layers*2];
        vector <double> structure;


        structure = var2struct(x);

        for (unsigned int i=0;i<structure.size();i++){
            array[i] = structure[i];
        }

        J_result = real_radiator_extern_(array, number_layers);

        delete[] array;

	    return J_result;
}

vector<double> var2struct(const vector<double> &x){
        const int n = x.size();
        int number_layers=C_NUMBER_LAYERS;
        int number_active_layers=C_NUMBER_ACTIVE_LAYERS;
        int layer_block,rest_layers,number_blocks;
        vector <double> array(2*number_layers,0);

        number_blocks=(number_active_layers+1);
        layer_block=(number_layers-number_active_layers)/number_blocks;
        rest_layers=number_layers-(layer_block*number_blocks+number_active_layers);
        if(rest_layers <0) {
            cout << "Error assigning values to array\n";
            exit(1);
        }

        int active_layer=0; // Active_layer and j have always the same values, but it is easier to follow the algorithm in this way
        int k=0;
        for (int j=0; j<number_blocks;j++){

        for (int i=j*layer_block+active_layer;i<(j+1)*layer_block+active_layer;i++){
            array[i] = x[k];
            k++;
        }
        if(active_layer < number_active_layers) array[(j+1)*layer_block+active_layer]=-1.0;
        active_layer++;
        }
        //Rest of layers
        if (rest_layers != 0) {
            for (int i=number_blocks*layer_block+number_active_layers; i<number_blocks*layer_block+number_active_layers+rest_layers; i++){
                array[i] = x[k];
                k++;

            }
        }

        for (int i=number_layers;i<2*number_layers;i++){
            array[i] = x[k];
            k++;
        }

        return array;
}
