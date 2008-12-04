
#include "rad_objfun.h"
#include <iostream>
#include <math.h>

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
