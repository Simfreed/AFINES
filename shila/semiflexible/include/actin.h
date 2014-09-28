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
#include "vector"

//=====================================
//actin filament class
class actin
{
    public:
        actin(double xcm, double ycm, double angle, double len, double fovx, double fovy, int nx, int ny, double vis);
        
        ~actin();
    
        double get_distance(double xp, double yp);

        double* get_intpoint(double xp, double yp);

        double get_int_angle(double xp, double yp);

        double* get_direction();

        double get_length();

        double get_xcm();
        
        double get_ycm();
        
        double get_angle();
        
        double* get_forces();

        void update_force(double f1, double f2, double f3);

        double* get_friction();
        
        double * get_start();
        
        double * get_end();

        std::vector<std::vector<int> > get_quadrants();
    
    private:
        double x,y,phi,ld, start[2], end[2], e[2], n[2], forces[3];
        
        double diameter, a_vis;
        
        std::vector<std::vector<int> > quad; //vector of two vectors(x and y quadrants) of integers
        
        std::vector<int> tmp;
};

#endif
