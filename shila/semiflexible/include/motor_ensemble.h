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

        motor_ensemble(double mdensity, double fovx, double fovy, double delta_t, double temp, double mlen, filament_ensemble_type* network, double v0,
                double stiffness, double ron, double roff, double rend, double actin_len, double vis, std::vector<double *> positions);

        ~motor_ensemble();

        void motor_walk(double t);

        void reshape();

        void motor_write(std::ofstream& fout);

        void motor_tension(std::ofstream& fout);

        void add_motor(motor<filament_ensemble_type> * m);

        void set_shear(double g);

    private:

        double fov[2], mrho, mld, mang, motorx, motory, alpha, gamma;
        int nm, s[2];
        filament_ensemble_type *f_network;
        std::vector<motor<filament_ensemble_type> *> n_motors;  
        std::string color;
};

#endif
