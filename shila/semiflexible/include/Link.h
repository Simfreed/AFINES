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

//=====================================
//Link class
class Link
{
    public:
        Link(double len, double stiffness, double bending_stiffness, filament* f, int aindex0, int aindex1);
        ~Link();

        double* get_hx();
        double* get_hy();
        double get_kb();
        double get_kl();
        double get_length();
        double get_posx();
        double get_posy();
        std::string get_color();
        std::string to_string();
        void step();
        void filament_update();
        bool operator==(const Link& that);    
        double get_stretch_force();

    protected:

        double hx[2],hy[2], phi, ld, stretch, forcex[2], forcey[2], torque[2], force_par[2],force_perp[2], kl, kb, xcm, ycm;
        
        int aindex[2];
        
        std::string color;
                
        filament *fil;
};

class MidLink : public Link
{
    public:
        void step();
};
#endif
