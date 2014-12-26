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

//=====================================
//actin rod class
class actin
{
    public:
        actin(double xcm, double ycm, double angle, double len, double fovx, double fovy, int nx, int ny, double vis); 
        ~actin();
    
        void update();

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

        void set_xcm(double xcm);

        void set_ycm(double ycm);

        void set_phi(double theta);

        void set_gay_berne(double sigma0, double eps0, double epsS, double epsE, double m, double n);
        
        double get_sigma0();
        
        double get_eps0();
        
        double get_chi();
        
        double get_chiPrime();

        double * calc_gay_berne(actin * a1);

        std::vector<std::vector<int> > get_quadrants();
       
        std::string write();
        
        std::string to_string();
        
        bool operator==(const actin& that);    
    

    private:
        double x,y,phi,ld, fov[2], nq[2], start[2], end[2], e[2], n[2], forces[3];
        
        double diameter, a_vis;
       
        double sigma0, eps0, mu, nu, chi, chiPrime;
        
        std::vector<std::vector<int> > quad; //vector of two vectors(x and y quadrants) of integers
        
        std::vector<int> tmp;
};

#endif
