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

#ifndef ANN_TB_NEURALNETWORK_H
#define ANN_TB_NEURALNETWORK_H

#include <vector>

namespace ann_toolbox {

class neural_network {
public:
	neural_network(unsigned int input_nodes_, unsigned int output_nodes_);
	virtual ~neural_network();

//	virtual void set_inputs(std::vector<double> &inputs);
	virtual void set_weights(const std::vector<double> &chromosome); 
    virtual const std::vector<double> compute_outputs(std::vector<double> &inputs) = 0;

	//virtual void SimulationStep(unsigned n_step_number, double f_time, double f_step_interval);

	unsigned int get_number_of_input_nodes() 	{ return m_inputs; };
	unsigned int get_number_of_output_nodes()	{ return m_outputs; };
	unsigned int get_number_of_weights() 		{ return m_weights.size(); };

protected:
	const char*		m_name;
	unsigned int	m_inputs, m_outputs;
	std::vector<double>	m_weights;
};

}
#endif
