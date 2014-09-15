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
class actin_ensemble;

//=====================================
//included dependences
#include "string"

//=====================================
//Link class
class Link
{
    public:
        Link(double len, double stiffness, actin_ensemble* network, int aindex0, int aindex1, std::string col);
        ~Link();

        double* get_heads();

        std::string get_color();

        void step();
        void actin_update();

    protected:

        double hx[2],hy[2], phi, ld, stretch, forcex[2], forcey[2], torque[2], force_par[2],force_perp[2], lk;
        
        int aindex[2];
        
        std::string color;
        
        actin_ensemble *actin_network;

};

class MidLink : public Link
{
    public:
        void step();
};

class BendingLink : public Link
{
    public:
        BendingLink(double len, double stiffness, double bending_stiffness, actin_ensemble* network, int aindex0, int aindex1, std::string col);

        void actin_update();

    protected:
        double bk;

};
#endif
