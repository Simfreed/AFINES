/*
 * motor.h
 *  
 *
 *  Created by Shiladitya Banerjee on 9/3/13.
 *  Copyright 2013 University of Chicago. All rights reserved.
 *
 */

//=====================================
//include guard
#ifndef __MOTOR_H_INCLUDED__
#define __MOTOR_H_INCLUDED__

//=====================================
// forward declared dependencies
//template <class filament_ensemble_type> class filament_ensemble;

//=====================================
//included dependences
#include "map"
#include "string"

//motor class
template <class filament_ensemble_type>
class motor 
{
    public:

        motor(double mx, double my, double mang, double mlen, filament_ensemble_type* network, int state0, int state1, 
                int findex0, int findex1, int rindex0, int rindex1, 
                double fovx, double fovy, double delta_t, double v0, double temp, double stiffness, double ron, double roff,
                double rend, double actin_len, double vis, string col);

        ~motor();

        string get_color();

        double tension();

        void attach(int hd);

        void brownian(double t, double gamma);

        void step_onehead(int hd);

        void step_twoheads();

        void actin_update();

        void update_shape();
        
        array<int,2> get_f_index();
        
        array<int,2> get_l_index();
        
        array<double,2> get_pos_a_end();
        
        array<int, 2> get_states();
        
        array<double, 2> get_hx();

        array<double, 2> get_hy();

        void detach_head(int hd);

        void move_end_detach(int hd, double pos);

        inline void reflect(double t, double gamma, double x1, double x2, double y1, double y2);

    private:

        double mphi,mld, mobility, onrate, stretch, vs, pos_temp, dm, fmax, mk, kon, koff, kend, dt, temperature;
        
        array<double,2> hx, hy, xm, ym, pos_a_end, fov;
        
        array<int,2> state, f_index, l_index;
        
        map<vector<int>, double> dist;
        
        string color;
        
        filament_ensemble_type *actin_network;
        
};

#endif
