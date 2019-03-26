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
#ifndef __MOTOR_ENSEMBLE_H_INCLUDED__
#define __MOTOR_ENSEMBLE_H_INCLUDED__

//=====================================
// forward declared dependencies
//class filament_ensemble;

//=====================================
//included dependences
#include "motor.h"
#include "spacer.h"

//motor ensemble class
template <class motor_type>
class motor_ensemble
{
    public:

        motor_ensemble(double mdensity, array<double, 2> myfov, double delta_t, double temp, double mlen, filament_ensemble* network, double v0,
                double stiffness, double max_ext_ratio, 
                double ron, double roff, double rend, 
                double fstall, double rcut,
                double vis, vector<array<double,3> > positions, string BC);

        motor_ensemble(vector<vector<double> > motors, array<double, 2> myfov, double delta_t, double temp, 
                double mlen, filament_ensemble * network, double v0, double stiffness, double max_ext_ratio, 
                double ron, double roff, double rend, 
                double fstall, double rcut,
                double vis, string BC); 
        
        motor_ensemble();
        
        ~motor_ensemble();

        int get_nmotors();

        void check_broken_filaments();

        void motor_walk(double t);

        void motor_update();
        
        void update_energies();
        
        double get_potential_energy();

        double get_kinetic_energy(); 

        void motor_write(ostream& fout);
        
        void motor_write_doubly_bound(ostream& fout);

        void print_ensemble_thermo();
        
        void motor_tension(ofstream& fout);

        void add_motor(motor_type * m);

        void set_shear(double g);
        
        void kill_heads(int i);
        
        void set_binding_two(double, double, double);

    protected:

        double mld, gamma, tMove;
        double ke, pe, v;

        array<double, 2> fov;
        filament_ensemble *f_network;
    
        vector<motor_type *> n_motors;  
};

class spacer_ensemble : public motor_ensemble<spacer>
{
    public:
        using motor_ensemble<spacer>::motor_ensemble;    
        void set_bending(double, double);

};

class xlink_ensemble : public motor_ensemble<motor>
{
    public:
        using motor_ensemble<motor>::motor_ensemble;    

};

#endif
