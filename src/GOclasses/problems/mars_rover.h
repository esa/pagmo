/****************************************************************
 * Copyright (C) 2009 Martin Peniak                             *
 * CopyPolicy: Released under the terms of the GNU GPL v2.0.    *
 ****************************************************************/

#ifndef _MARSROVER_H
#define _MARSROVER_H

#include <string>
#include <vector>
#include "../problems/GOproblem.h"
#include "../basic/population.h"
using namespace std;

class __PAGMO_VISIBLE mars_rover: public GOProblem
{
        public:
	///Constructor
        mars_rover(int dim, double bounds):GOProblem(vector<double>(dim, -bounds), vector<double>(dim, bounds)),num_hidden(0),recurr(0),num_steps(5000),num_trials(10),bin_sensors(1),graphics(0),max_velocity(3.0),rand_pos_range(23),kill_threshold(0.0),step_size(0.1),full_screen(0),window_width(1280),window_height(960),pos_x(0.0),pos_y(1.2),pos_z(0.0),terrain_scale(1),surface_slip1(0.0008),surface_slip2(0.0008),erp(0.2),cfm(0.00015),island_terrain(8,0),island_faults(8,0),infrared_sensors(1),active_vision(1),proprioception(0),num_pixels(5),rgb(1),zoomable(0),max_pan(45),max_tilt(22.5),pan_speed(1),tilt_speed(1),debug(0) {}
       
       ///Setters
       void set_num_hidden(int num_hidden_) {
	num_hidden = num_hidden_;
       }

        virtual void pre_evolution(Population &) const;

        // Cloning method
        virtual mars_rover *clone() const;

        // Object descriptor
        virtual string id_object() const;

        private:
        //Objective function (i.e. fitness function)
	int    num_hidden;
        int    recurr;
        int    num_steps;
        int    num_trials;
        int    bin_sensors;
        int    graphics;
        int    full_screen;
        int    window_width;
        int    window_height;
        int    infrared_sensors;
        int    active_vision;
        int    proprioception;
        int    num_pixels;
        int    rgb;
        int    zoomable;
        int    debug;
        float  max_velocity;
        float  rand_pos_range;
        float  kill_threshold;
        float  step_size;
        float  pos_x;
        float  pos_y;
        float  pos_z;
        float  terrain_scale;
        float  max_pan;
        float  max_tilt;
        float  pan_speed;
        float  tilt_speed;
        double surface_slip1;
        double surface_slip2;
        double erp;
        double cfm;
        vector<int> island_terrain;
        vector<int> island_faults;
        virtual double objfun_(const std::vector<double> &x, size_t m_id) const;

        mutable size_t          m_random_seed;
};

#endif
