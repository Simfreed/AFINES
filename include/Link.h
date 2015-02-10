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

        double* get_hx();
        
        double* get_hy();
        
        double get_kb();
        
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
        
        double * get_forces();

    protected:

        double hx[2],hy[2], phi, ld, stretch, forcex[2], forcey[2], torque[2], force_par[2],force_perp[2], kl, kb, xcm, ycm;
        
        int aindex[2];
        
        filament *fil;
};

class MidLink : public Link
{
    public:
        MidLink(double len, double stiffness, double bending_stiffness, filament* f, int aindex0, int aindex1);
        
        void step();
        
        double get_stretch_force();
        
        void set_linear(bool lin);

    protected:

        bool is_linear;
};
#endif
