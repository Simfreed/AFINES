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
#include "globals.h"

//motor class
template <class filament_ensemble_type>
class motor 
{
    public:

        motor(array<double, 3> pos, double mlen, filament_ensemble_type* network, 
                array<int, 2> mystate, array<int, 2> myfindex, array<int, 2> myrindex,
                array<double, 2> myfov, double delta_t, double v0, double temp, double stiffness, double ron, double roff,
                double rend, double actin_len, double vis, string BC);

        ~motor();

        string get_BC();

        double tension();

        void attach(double t, int hd);

        void update_force(double t);
        
        void update_position(double t);

        array<double, 2> boundary_check(int i, double t, double vx, double vy);

        void step_onehead(int hd);

        void step_twoheads();

        void actin_update();

        void update_shape();
        
        void set_shear(double g);
        
        array<int,2> get_f_index();
        
        array<int,2> get_l_index();
        
        array<double,2> get_pos_a_end();
        
        array<int, 2> get_states();
        
        array<double, 2> get_hx();

        array<double, 2> get_hy();

        void detach_head(int hd);

        void move_end_detach(int hd, double pos);

        inline void reflect(double t, double gamma, double x1, double x2, double y1, double y2);
        
        inline void periodic(double t, double gamma, double x1, double x2, double y1, double y2);

    private:

        double mphi,mld, mobility, vs, dm, fmax, mk, kon, koff, kend, dt, temperature, 
               actin_damp, force, damp, shear;
        
        array<double,2> hx, hy, xm, ym, pos_a_end, fov, prv_rnd_x, prv_rnd_y;
        
        array<int,2> state, f_index, l_index;
        
        map<vector<int>, double> dist;
        
        string BC;
        
        filament_ensemble_type* actin_network;
        
};

#endif
