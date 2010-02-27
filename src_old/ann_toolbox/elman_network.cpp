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

#include "elman_network.h"

using namespace ann_toolbox;

// Constructor
elman_network::elman_network(unsigned int input_nodes_, unsigned int hidden_nodes_,
	unsigned int output_nodes_, const std::vector<double> &w) : 
    	neural_network(input_nodes_, output_nodes_),
		m_hidden(hidden_nodes_)
{
	// the number of weights is equal to all the inputs (and a bias)
	// for every hidden node, plus the connections from every hidden
	// node to every output, i.e. it is fully connected. The feedback
	// nodes are connected to the hidden (middle) layer too.
	unsigned int wghts = (m_inputs + 1 + m_hidden) * m_hidden + (m_hidden + 1) * m_outputs;
	m_weights = std::vector<double>(wghts, 0);
	
	// the memory (feedback) values are stored in this vector
	m_memory  = std::vector<double>(m_hidden, 0);
	
	if(! w.empty()) set_weights(w);

}

// Destructor
elman_network::~elman_network() {}

// Computing the outputs
const std::vector<double> elman_network::compute_outputs(std::vector<double> &inputs) 
{
	// check for correct input size
	if(inputs.size() != m_inputs) {
		pagmo_throw(value_error, "incorrect size of input vector");
	}
	
	// generate values for the hidden nodes
	std::vector<double> hidden(m_hidden, 0);
	unsigned int i = 0, j, offset;
    for(  ; i < m_hidden; i++ ) {
		offset = (m_inputs + m_hidden + 1) * i;
        // Set the bias
        hidden[i] = m_weights[offset + m_inputs];
        
        for( j = 0; j < m_inputs; j++ ) {
            // Add the weighted input
            hidden[i] += m_weights[offset + j] * inputs[j];			 
        }

        for( j = 0; j < m_hidden; j++ ) {
            // Add the weighted memory
            hidden[i] += m_weights[offset + m_inputs + 1 + j] * m_memory[j];			 
        }
        
        // Apply the transfer function (a sigmoid with output in [0,1])
        hidden[i] = 1.0/( 1 + exp( -hidden[i] ));	  
    }

    // Copy the hidden to the memory:
	m_memory = hidden;

	// generate values for the output nodes	
	std::vector<double> outputs(m_outputs, 0);
    for( i = 0; i < m_outputs; i++ ) {
    	// Offset for the weights to the output nodes
 		offset = (m_inputs + m_hidden + 1) * m_hidden + (m_hidden + 1) * i;

        // add the bias (wheigted by the first wheigt to the i^th output node
        outputs[i] = m_weights[offset + m_hidden];		 

		for ( j = 0; j < m_hidden; j++) {
			outputs[i] += m_weights[offset + j] * hidden[j];
	    }
	    
		// Apply the transfer function (a sigmoid with output in [0,1])
		outputs[i] = 1.0/(1 + exp( -outputs[i]));
    }

    return outputs;
}