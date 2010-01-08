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

#ifndef ANN_TB_CTRNN_H
#define ANN_TB_CTRNN_H

#include "neural_network.h"

namespace ann_toolbox {
	
/**
 * A Continuous-Time Recurrent Neural Networks (CTRNNs) is a dynamic neural network
 * allowing recurrent connections and considers continuous time as a factor.
 * TODO
 * More info: add link
 */	
class ctrnn : public neural_network {
public:
	/// Constructor
	/**
	 * Creates a new ctrnn (continuous-time recurrent neural network) object, which is derived
	 * from the neural_network class.
	 * \param input_nodes	the number of input nodes
	 * \param hidden_nodes	the number of nodes in the hidden layer
	 * \param output_nodes	the number of output nodes (default = 1)
	 * \param w				the weights, with which the neural network is initiated (empty by default)
	 * \return a perceptron object
	 */
	ctrnn(unsigned int input_nodes_, unsigned int hidden_nodes_, 
			unsigned int output_nodes_ = 1, const std::vector<double> &w = std::vector<double>());	

	/// Destructor
    ~ctrnn();

	/// Getter/SetterFunctions
	virtual void set_weights(const std::vector<double> &w); 
	//void set_time_step(double ts);

	/// Compute Outputs
	const std::vector<double> compute_outputs(std::vector<double> &inputs);

protected:
	// number of hidden nodes
	unsigned int	m_hidden;
	double	m_time_step;
	// a vector to store the memory of the network (feedback nodes)
	std::vector<double>	m_hidden_neurons;
	std::vector<double>	m_output_neurons;	
};

}
#endif

