/*
 *  actin.cpp
 *  
 *
 *  Created by Shiladitya Banerjee on 9/3/13.
 *  Copyright 2013 University of Chicago. All rights reserved.
 *
 */

//=====================================
//include guard
#ifndef __ACTIN_H_INCLUDED__
#define __ACTIN_H_INCLUDED__

//=====================================
// forward declared dependencies

//=====================================
//included dependences
#include "string"
#include "vector"
#include "globals.h"

//=====================================
//actin rod class
class actin
{
    public:

        actin();
        
        actin(double xcm, double ycm, double len, double vis); 
        
        actin(double xcm, double ycm, double vx, double vy, double len, double vis); 
        
        actin(const actin& other);
        
        ~actin();
    
        void update();

        double get_length();

        double get_xcm();
        
        double get_ycm();
        
        void update_velocity(double vx, double vy);
        
        void update_force(double f1, double f2);
        
        void reset_velocity();
        
        void reset_force();

        array<double,2> get_velocity();
        
        array<double,2> get_force();

        double get_friction();
        
        double get_viscosity();

        double get_vsquared();
        
        void set_xcm(double xcm);

        void set_ycm(double ycm);

        string write();
        
        string to_string();
        
        bool operator==(const actin& that);    
        

    private:
        
        double x, y, ld, a_vis, friction;

        array<double, 2> force, velocity;
        
};

#endif
