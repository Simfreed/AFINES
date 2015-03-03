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
        
        actin(double xcm, double ycm, double len, double fovx, double fovy, int nx, int ny, double vis); 
        
        actin(const actin& other);
        
        ~actin();
    
        void update();

        double get_distance(double xp, double yp);


        double get_int_angle(double xp, double yp);


        double get_length();

        double get_xcm();
        
        double get_ycm();
        
        double get_angle();
        
        void update_force(double f1, double f2, double f3);

        array<double,2> get_fov();

        array<int,2> get_nq();
        
        array<double,2> get_intpoint(double xp, double yp);
        
        array<double,2> get_forces();

        double get_friction();
        
        double get_viscosity();

        void set_xcm(double xcm);

        void set_ycm(double ycm);

        vector<vector<int> > get_quadrants();
       
        string write();
        
        string to_string();
        
        bool operator==(const actin& that);    
        

    private:
        
        double x, y, ld, a_vis, friction;

        array<double,2> fov; 
        
        array<int,2> nq;
        
        array<double, 2> forces;
        
        vector<vector<int> > quad; //vector of two vectors(x and y quadrants) of integers
        
        vector<int> tmp;
};

#endif
