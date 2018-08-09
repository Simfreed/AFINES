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
class filament_ensemble;

//=====================================
//included dependences
#include "globals.h"

//motor class
class motor 
{
    public:

        motor(array<double, 3> pos, double mlen, filament_ensemble* network, 
                array<int, 2> mystate, array<int, 2> myfindex, array<int, 2> myrindex,
                array<double, 2> myfov, double delta_t, double v0, double temp, double stiffness, double max_ext_ratio, 
                double ron, double roff, double rend, 
                double fstall, double rcut,
                double vis, string BC);
        
        motor(array<double, 4> pos, double mlen, filament_ensemble* network, 
                array<int, 2> mystate, array<int, 2> myfindex, array<int, 2> myrindex,
                array<double, 2> myfov, double delta_t, double v0, double temp, double stiffness, double max_ext_ratio,
                double ron, double roff, double rend, 
                double fstall, double rcut,
                double vis, string BC);
        
        motor();
        
        virtual ~motor();

        string get_BC();

        virtual bool allowed_bind( int hd, array<int, 2> fl_idx);
        
        bool attach( int hd);

        void relax_head( int hd);

        void kill_head( int hd);
        
        void identify();

        void update_position_attached(int hd);
        
        virtual void update_force();
        
        virtual void update_force_fraenkel_fene();
        
        void update_angle();
        
        virtual void brownian_relax(int hd);

        array<double, 2> boundary_check(int i,  double vx, double vy);

        void step_onehead( int hd);

        void step_twoheads();

        void filament_update_hd(int hd, array<double, 2> f);
        
        void filament_update();

        void update_shape();
        
        void set_shear(double g);
        
        array<int,2> get_f_index();
        
        array<int,2> get_l_index();
        
        void init_l_index(int hd, int idx);
        
        void set_l_index(int hd, int idx);
        
        void set_f_index(int hd, int idx);
        
        array<double,2> get_pos_a_end();
        
        void set_pos_a_end(int hd, double pos);
        
        double get_pos_a_end(int hd);
        
        array<double,2> get_force();
        
        array<int, 2> get_states();
        
        array<double, 2> get_hx();

        array<double, 2> get_hy();

        void detach_head(int hd, array<double, 2> pos);

        void detach_head_without_moving(int hd);
        
        void update_pos_a_end(int hd, double pos);

        inline void reflect(double t, double gamma, double x1, double x2, double y1, double y2);
        
        inline void periodic(double t, double gamma, double x1, double x2, double y1, double y2);

        double get_stretching_energy();
        
        double get_stretching_energy_fene();

        double get_kinetic_energy();

        double get_angle();
        
        virtual double metropolis_prob(int hd, array<int, 2> fl_idx, array<double, 2> newpos, double maxprob);

        array<double, 2> generate_off_pos(int hd);

        string to_string();
        
        string write();
        
        void remove_from_spring(int hd);

        void add_to_spring(int hd);
        
        void inc_l_index(int hd);
    
        void set_binding_two(double ron2, double roff2, double rend2);

    public:

        double mld, vs, stall_force, max_bind_dist, max_bind_dist_sq, mk, kon, koff, kend, dt, temperature, 
               damp, shear, max_ext, eps_ext, kinetic_energy, bd_prefactor, tension, len,
               kon2, koff2, kend2;
        
        array<double,2> hx, hy, pos_a_end, fov, prv_rnd_x, prv_rnd_y, force, disp, direc;

        array<array<double, 2>, 2> ldir_bind, bind_disp;

        array<int,2> state, f_index, l_index, spring_mot_idx;
        
        map<vector<int>, double> dist;
        array<bool, 2> at_barbed_end;
        
        string BC;
        
        filament_ensemble* filament_network;
        
};

#endif
