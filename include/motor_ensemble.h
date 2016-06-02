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
#include "vector"

//motor ensemble class
template <class filament_ensemble_type>
class motor_ensemble
{
    public:

        motor_ensemble(double mdensity, array<double, 2> myfov, double delta_t, double temp, double mlen, filament_ensemble_type* network, double v0,
                double stiffness, double max_ext_ratio, 
                double ron, double roff, double rend, 
                double fstall, double fbreak, double engBind,
                double vis, vector<array<double,3> > positions, string BC);

        motor_ensemble(vector<vector<double> > motors, array<double, 2> myfov, double delta_t, double temp, 
                double mlen, filament_ensemble_type * network, double v0, double stiffness, double max_ext_ratio, 
                double ron, double roff, double rend, 
                double fstall, double fbreak, double engBind,
                double vis, string BC); 
        
        ~motor_ensemble();

        int get_nmotors();

        void check_broken_filaments();

        void motor_walk(double t);

        void motor_update();
        
        void update_energies();
        
        double get_potential_energy();

        void motor_write(ostream& fout);

        void print_ensemble_thermo();
        
        void motor_tension(ofstream& fout);

        void add_motor(motor<filament_ensemble_type> * m);

        void set_shear(double g);
        
        void kill_heads(int i);

    private:

        double mld, gamma, tMove;
        double ke, pe, v;

        array<double, 2> fov;
        filament_ensemble_type *f_network;
        vector<motor<filament_ensemble_type> *> n_motors;  
};

#endif
