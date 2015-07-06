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
                double stiffness, double ron, double roff, double rend, double actin_len, double vis, vector<array<double,3> > positions, string BC);

        ~motor_ensemble();

        int get_nmotors();

        void check_broken_filaments();

        void motor_walk(double t);

        void reshape();

        void motor_write(ostream& fout);

        void motor_tension(ofstream& fout);

        void add_motor(motor<filament_ensemble_type> * m);

        void set_shear(double g);

    private:

        double mrho, mld, mang, motorx, motory, alpha, gamma, tMove;
        int nm;
        
        array<double, 2> fov;
        array<int, 2> s;
        filament_ensemble_type *f_network;
        vector<motor<filament_ensemble_type> *> n_motors;  
        string color;
};

#endif
