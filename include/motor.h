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
                array<double, 2> myfov, double delta_t, double v0, double temp, double stiffness, double max_ext_ratio, 
                double ron, double roff,
                double rend, double actin_len, double vis, string BC);
        
        motor(array<double, 4> pos, double mlen, filament_ensemble_type* network, 
                array<int, 2> mystate, array<int, 2> myfindex, array<int, 2> myrindex,
                array<double, 2> myfov, double delta_t, double v0, double temp, double stiffness, double max_ext_ratio,
                double ron, double roff,
                double rend, double actin_len, double vis, string BC);

        ~motor();

        string get_BC();

        double tension();

        bool attach( int hd);

        void relax_head( int hd);

        void kill_head( int hd);
        
        void update_force();
        
        void update_force_fraenkel_fene();
        
        void update_angle();
        
        void brownian_relax(int hd);

        array<double, 2> boundary_check(int i,  double vx, double vy);

        void step_onehead( int hd);

        void step_twoheads();

        void actin_update();

        void update_shape();
        
        void set_shear(double g);
        
        array<int,2> get_f_index();
        
        array<int,2> get_l_index();
        
        array<double,2> get_pos_a_end();
        
        array<double,2> get_force();
        
        array<int, 2> get_states();
        
        array<double, 2> get_hx();

        array<double, 2> get_hy();

        void detach_head(int hd);

        void move_end_detach( int hd, double pos);

        inline void reflect(double t, double gamma, double x1, double x2, double y1, double y2);
        
        inline void periodic(double t, double gamma, double x1, double x2, double y1, double y2);

        double get_stretching_energy();
        
        double get_stretching_energy_fene();

        double get_kinetic_energy();
        
        string to_string();
        
        string write();
    private:

        double mphi,mld, mobility, vs, dm, stall_force, mk, kon, koff, kend, dt, temperature, 
               actin_damp, damp, shear, max_ext, eps_ext, kinetic_energy, bd_prefactor;
        
        array<double,2> hx, hy, xm, ym, pos_a_end, fov, prv_rnd_x, prv_rnd_y, force, disp;
        
        array<int,2> state, f_index, l_index;
        
        map<vector<int>, double> dist;
        
        string BC;
        
        filament_ensemble_type* actin_network;
        
};

#endif
