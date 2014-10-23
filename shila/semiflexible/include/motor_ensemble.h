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
//class actin_ensemble;

//=====================================
//included dependences
#include "motor.h"
#include "vector"

//motor ensemble class
class motor_ensemble
{
    public:

        motor_ensemble(double mdensity, double fovx, double fovy, double delta_t, double mlen, actin_ensemble* network, double v0,
                double stiffness, double ron, double roff, double rend, double actin_len, double vis, std::vector<double *> positions);

        ~motor_ensemble();

        void motor_walk();

        void reshape();

        void motor_write(std::ofstream& fout);

        void motor_tension(std::ofstream& fout);

        void add_motor(motor * m);

    private:

        double fov[2], mrho, mld, mang, motorx, motory, alpha;
        int nm, s[2];
        actin_ensemble *a_network;
        std::vector<motor *> n_motors;  
        std::string color;
};

#endif
