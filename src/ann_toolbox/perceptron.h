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

#ifndef ANN_TB_PERCEPTRON_H
#define ANN_TB_PERCEPTRON_H

#include "neural_network.h"

namespace ann_toolbox {
	
/**
 * A simple perceptron (a type of artificial neural network), representing the 
 * simplest kind of feedforward neural network. This basically refers to a linear
 * classifier. 
 * More info: http://en.wikipedia.org/wiki/Perceptron
 */	
class perceptron : public neural_network {
public:
	/// Constructors
	/**
	 * Creates a new perceptron object, which is derived from the neural_network
	 * class. If initial weights are handed over it calls the set_weights function
	 * to initalize the weights of the neural network.
	 * \param input_nodes	the number of input nodes
	 * \param output_nodes	the number of output nodes (not mandatory, default = 1)
	 * \param w				the weights, with which the neural network is initiated (not mandatory)
	 * \return a perceptron object
	 */
	perceptron(unsigned int input_nodes_, unsigned int output_nodes_ = 1, const std::vector<double> &w = std::vector<double>());

	/// Destructor
    ~perceptron();

	/// Compute Outputs
	const std::vector<double> compute_outputs(std::vector<double> &inputs);

protected:
};

}
#endif
