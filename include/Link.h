/*
 *  link.h
 *
 *  Created by Simon Freedman on 9/10/2014
 *  Copyright 2014 University of Chicago. All rights reserved.
 *
 */

//=====================================
//include guard
#ifndef __LINK_H_INCLUDED__
#define __LINK_H_INCLUDED__


//=====================================
// forward declared dependencies
class filament;

//=====================================
//included dependences
#include "string"
#include "vector"
#include "globals.h"

//=====================================
//Link class
class Link
{
    public:
        Link();
        
        Link(double len, double stiffness, double bending_stiffness, filament* f, int aindex0, int aindex1);
        
        virtual ~Link();

        array<double, 2> get_hx();
        
        array<double, 2> get_hy();
        
        double get_kl();
        
        double get_length();
        
        double get_xcm();
        
        double get_ycm();
        
        string to_string();
        
        string write();
        
        void step();
        
        void filament_update();
        
        bool operator==(const Link& that);    
        
        bool is_similar(const Link& that);    

        virtual double get_stretch_force();

        void set_aindex1(int i);
        
        array<double,4> get_forces();
        
        double get_distance(double xp, double yp);

        double get_int_angle(double xp, double yp);
        
        array<double,2> get_intpoint(double xp, double yp);
        
        vector<array<int,2> > get_quadrants();
       
        void quad_update();
        
        array<double, 2> get_direction();

    protected:

        double xcm, ycm, phi, l0, kl;
       
        array<double,2> hx, hy;
        
        array<int, 2> aindex;
        
        filament *fil;
        
        vector< array<int,2> > quad; //vector of two vectors(x and y quadrants) of integers
};
#endif
