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

// Created by Juxi Leitner on 2009-12-21.
// based on the TwoDee Artificial Neural Network Code

#include <cstdlib>
#include <cmath>
#include <vector>        
#include <exception>
#include "../exceptions.h"

#include "ctrnn.h"

using namespace ann_toolbox;

// create these funtions to easier access the weights
#define input_to_hidden_weights(idx)	m_weights[0 + (idx)]
#define hidden_to_hidden_weights(idx)	m_weights[m_inputs * m_hidden + (idx)]
#define hidden_bias(idx)				m_weights[m_inputs * m_hidden + m_hidden * m_hidden	+ (idx)]
#define hidden_taus(idx)				m_weights[m_inputs * m_hidden + m_hidden * m_hidden	+ m_hidden + (idx)]
#define hidden_to_output_weights(idx)	m_weights[m_inputs * m_hidden + m_hidden * m_hidden	+ m_hidden + m_hidden + (idx)]
#define output_bias(idx)				m_weights[m_inputs * m_hidden + m_hidden * m_hidden	+ m_hidden + m_hidden + m_hidden * m_outputs + (idx)]

// Constructor
ctrnn::ctrnn(unsigned int input_nodes_, unsigned int hidden_nodes_, unsigned int output_nodes_,
	const std::vector<double> &w) : 
    	neural_network(input_nodes_, output_nodes_),
		m_hidden(hidden_nodes_),
		m_time_step(0.1),
		m_hidden_neurons(std::vector<double>(hidden_nodes_, 0)),
		m_output_neurons(std::vector<double>(output_nodes_, 0))
{
	// the number of weights
	unsigned int wghts = m_inputs * m_hidden// synaptic weights from input to hidden layer
					+ m_hidden * m_hidden	// synaptic weights from hidden to hidden layer
					+ m_hidden + m_hidden	// bias and taus of the hidden layer
					+ m_hidden * m_outputs	// synaptic weights from hidden to output layer
					+ m_outputs;			// bias of the output layer

	m_weights = std::vector<double>(wghts, 0);
	
	lowbound_weights = -10.0;
    uppbound_weights =  10.0;
    lowbound_bias    = -10.0;
	uppbound_bias    =  10.0;
    lowbound_tau     = -1.0;
    uppbound_tau     =  2.0;	
	
	if(! w.empty()) set_weights(w);
}

// Destructor
ctrnn::~ctrnn() {}

void ctrnn::set_weights(const std::vector<double> &w) {
	if(w.size() != m_weights.size()) {
		pagmo_throw(value_error, "number of weights is incorrect");
	}
	m_weights = w;
	    
	// scale them according to the boundary values (maybe need a set funciton for those?)
	unsigned int idx = 0, i;
	
	// Scale the weigth values to be in the right range
    for (i = 0; i < m_inputs * m_hidden; i++) {
        input_to_hidden_weights(i) = w[idx++] * (uppbound_weights - lowbound_weights) + lowbound_weights;
    }

    for (i = 0; i < m_hidden * m_hidden; i++) {
        hidden_to_hidden_weights(i) = w[idx++] * (uppbound_weights - lowbound_weights) + lowbound_weights;
    }

    for (i = 0; i < m_hidden; i++) {
        hidden_bias(i) = w[idx++] * (uppbound_bias - lowbound_bias) + lowbound_bias;
    }

    for (i = 0; i < m_hidden; i++) {
        hidden_taus(i) = pow(10, (lowbound_tau + ((uppbound_tau - lowbound_tau) * w[idx++])));
    }

    for (i = 0; i < m_hidden * m_outputs; i++) {
        hidden_to_output_weights(i) = w[idx++] * (uppbound_weights - lowbound_weights) + lowbound_weights;
    }

    for (i = 0; i < m_outputs; i++) {
        output_bias(i) = w[idx++] * (uppbound_bias - lowbound_bias) + lowbound_bias;
    }

	// set hidden and output layers 0
	std::fill(m_hidden_neurons.begin(), m_hidden_neurons.end(), 0.0);
	std::fill(m_output_neurons.begin(), m_output_neurons.end(), 0.0);	
	
}

// Computing the outputs
const std::vector<double> ctrnn::compute_outputs(std::vector<double> &inputs) {
	// check for correct input size
	if(inputs.size() != m_inputs) {
		pagmo_throw(value_error, "incorrect size of input vector");
	}
	
	unsigned int i, j;
	// Update delta state of hidden layer from inputs:
	std::vector<double> tempH(m_hidden, 0);
	for(i = 0; i < m_hidden; i++) {
    	tempH[i] = - m_hidden_neurons[i];
  
    	for( j = 0; j < m_inputs; j++) {
        	// weight * sigmoid(state)
			tempH[i] += input_to_hidden_weights(i * m_inputs + j) * inputs[j] ;	  
		}	  
	}

	double h;
	// Update delta state from hidden layer, self-recurrent connections:
	for(i = 0; i < m_hidden; i++) {
    	for(j = 0; j < m_hidden; j++) {
        	h = (double(1.0)/( exp(-( m_hidden_neurons[j] + hidden_bias(j))) + 1.0 ));
        	tempH[i] += hidden_to_hidden_weights(i * m_hidden + j) * h;
    	}
	}
	
	// Update the hidden layer
	for( i = 0; i < m_hidden; i++) {
    	m_hidden_neurons[i] += tempH[i] * m_time_step/hidden_taus(i);
	}

	// Update the output layer
	for(i = 0; i < m_outputs; i++) {
    	for(j = 0; j < m_hidden; j++) {
        	h = (double(1.0)/( exp(-( m_hidden_neurons[j] + hidden_bias(j))) + 1.0 ));
        	m_output_neurons[i] += hidden_to_output_weights(i * m_hidden + j) * h;
    	}

    	// Compute the activation function immediately, since this is
    	// what we return and since the output layer is not recurrent:
    	m_output_neurons[i] = double(1.0)/( exp(-( m_output_neurons[i] + output_bias(i))) + 1.0 );
	}	
	
    return m_output_neurons;
}