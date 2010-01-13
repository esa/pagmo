/****************************************************************
 * Copyright (C) 2009 Martin Peniak                             *
 * CopyPolicy: Released under the terms of the GNU GPL v2.0.    *
 ****************************************************************/

#include "../problems/mars_rover.h"
#include <boost/lexical_cast.hpp>
#include <cstdio>
#include <exception>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

//Objective function (i.e. fitness function); passed vector containing genes and island ID
double mars_rover::objfun_(const std::vector<double> &x, size_t m_id) const
{
    //Corrects Island ID to Island number (e.g. the 4th island will have the m_id of 9);
    m_id = (m_id-1)/2;
    cout<<"running island: "<<m_id<<" terrain ID: "<<island_terrain[m_id]<<" number of faults: "<<island_faults[m_id]<<endl;
    //Local variable for returning fitness
    double fitness = 0;

    //Holds genotype size
    const size_t size = x.size();

    //Variables used for passing genes as argumnets
    char buf[200];
    char parameters[100000];

    //Output file name passed to the simulator to save the fitness
    char *output_name;
    output_name = tempnam(0,boost::lexical_cast<string>(m_random_seed).c_str());

    //Create a string with all the parameters
    strcpy(parameters,"./ESA_Project ");
    sprintf(buf,"%i ",num_hidden);
    strcat(parameters, buf);
    sprintf(buf,"%i ",recurr);
    strcat(parameters, buf);
    sprintf(buf,"%i ",num_steps);
    strcat(parameters, buf);
    sprintf(buf,"%i ",num_trials);
    strcat(parameters, buf);
    sprintf(buf,"%i ",bin_sensors);
    strcat(parameters, buf);
    sprintf(buf,"%f ",max_velocity);
    strcat(parameters, buf);
    sprintf(buf,"%f ",rand_pos_range);
    strcat(parameters, buf);
    sprintf(buf,"%f ",kill_threshold);
    strcat(parameters, buf);
    sprintf(buf,"%f ",step_size);
    strcat(parameters, buf);

    //Add the random seed
    sprintf(buf,"%i ",m_random_seed);
    strcat(parameters,buf);

    //Add file path to the parameters
    sprintf(buf,"%s ",output_name);
    strcat(parameters,buf);

    //Graphics
    sprintf(buf,"%i ",graphics);
    strcat(parameters,buf);

    //Path to the terrain
    sprintf(buf,"%s%i ","Data/models/terrains/",island_terrain[m_id]);
    strcat(parameters,buf);

    //Number of failed sensors
    sprintf(buf,"%i ",island_faults[m_id]);
    strcat(parameters,buf);

    //Initial rover position, pos_y is its starting vertical position used throughout the evolution
    sprintf(buf,"%f ",pos_x);
    strcat(parameters,buf);
    sprintf(buf,"%f ",pos_y);
    strcat(parameters,buf);
    sprintf(buf,"%f ",pos_z);
    strcat(parameters,buf);

    //Terrain scale
    sprintf(buf,"%f ",terrain_scale);
    strcat(parameters,buf);

    //Surface x,y slip
    sprintf(buf,"%f ",surface_slip1);
    strcat(parameters,buf);
    sprintf(buf,"%f ",surface_slip2);
    strcat(parameters,buf);

    //Error Reduction Parameter and Constraint Force Mixing (ODE)
    sprintf(buf,"%f ",erp);
    strcat(parameters,buf);
    sprintf(buf,"%f ",cfm);
    strcat(parameters,buf);

    //Width and height of the window
    sprintf(buf,"%i ",window_width);
    strcat(parameters,buf);
    sprintf(buf,"%i ",window_height);
    strcat(parameters,buf);

    //Fullscreen mode flag
    sprintf(buf,"%i ",full_screen);
    strcat(parameters,buf);

    //Sensory inputs
    sprintf(buf,"%i ",infrared_sensors);
    strcat(parameters,buf);
    sprintf(buf,"%i ",active_vision);
    strcat(parameters,buf);
    sprintf(buf,"%i ",proprioception);
    strcat(parameters,buf);

    //Active vision details
    sprintf(buf,"%i ", num_pixels);
    strcat(parameters,buf);
    sprintf(buf,"%i ",rgb);
    strcat(parameters,buf);
    sprintf(buf,"%i ",zoomable);
    strcat(parameters,buf);
    sprintf(buf,"%f ",max_pan);
    strcat(parameters,buf);
    sprintf(buf,"%f ",max_tilt);
    strcat(parameters,buf);
    sprintf(buf,"%f ",pan_speed);
    strcat(parameters,buf);
    sprintf(buf,"%f ",tilt_speed);
    strcat(parameters,buf);

    //Mode
    sprintf(buf,"%i ",debug);
    strcat(parameters,buf);

    //Genotype size
    sprintf(buf,"%i ",size);
    strcat(parameters,buf);

    //Add genes to the parameters
    for (size_t i = 0; i < size; ++i)
    {
        sprintf(buf,"%f ",x.at(i));
        strcat(parameters,buf);
    }
    int status = system(parameters);
    ifstream res_stream(output_name);
    if (res_stream.is_open())
    {
        string line;
        getline(res_stream,line);
        fitness = atof(line.c_str());
    }
    else
    {
        status = 1;
    }

    // Try removing the temporary files.
    remove(output_name);
    if (status)
    {
        cout<<"ERROR: file containing fitness does not exist..."<<endl;
        cout<<"       returning zero fitness"<<endl;
    }

    //Returns negative fitness, due to PaGMO's minimisation
    return -fitness;
}

//Re-evaluate the population with respect to the new seed (Internal Sampling Method)
void mars_rover::pre_evolution(Population &pop) const
{
        size_t m_id = pop.m_id;
        m_random_seed = static_rng_uint32()();
        for (size_t i=0; i<pop.size(); ++i)
        {
          pop[i] = Individual(*this, pop[i].getDecisionVector(), m_id);
        }
}

//Cloning method: this will return a pointer to a new instance of the problem, copied from the current one
mars_rover *mars_rover::clone() const
{
        return new mars_rover(*this);
}

// Returns a description of the problem
std::string mars_rover::id_object() const
{
        return std::string("Mars rover simulation problem.");
}
