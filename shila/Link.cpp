/*
 *  Link.cpp
 *  
 *
 *  Created by Simon Freedman on 9/10/2014
 *  Copyright 2014 University of Chicago. All rights reserved.
 *
 */

#include "generic_functions.cpp"
#include "Link.h"
#include "actin_ensemble.h"

Link::Link(double len, double stiffness, actin_ensemble* network, 
        int aindex0, int aindex1, std::string col)
{
    lk              =   stiffness;
    ld              =   len;
    aindex[0]       =   aindex0;
    aindex[1]       =   aindex1;
    actin_network   =   network;

    // Set the coordinates of the heads:
    this->step();
    color           =   col; 

}


double* Link::get_heads()
{
    double h[4];
    double *gh;
    h[0]=hx[0];
    h[1]=hy[0];
    h[2]=hx[1];
    h[3]=hy[1];
    gh=h;
    return gh;
}

std::string Link::get_color()
{
    return color;
}

// stepping kinetics
void Link::step()
{

    //CONVENTION: head 0 will be connected to the POINTY end of a filament
    //            head 1 will be connected to the BARBED end of a filament
    hx[0]=actin_network->get_ends(aindex[0])[2];
    hy[0]=actin_network->get_ends(aindex[0])[3];
    hx[1]=actin_network->get_ends(aindex[0])[0];
    hy[1]=actin_network->get_ends(aindex[0])[1];

    phi=atan2(hy[1]-hy[0],hx[1]-hx[0]);

}

void Link::actin_update()
{
    stretch         =   dis_points(hx[0],hy[0],hx[1],hy[1])-ld;
    std::cout<<"DEBUG: actin_update: color = "<< color<< "\tstretch = "<<stretch<<"\n";

    forcex[0]       =   lk * stretch * cos(phi); //-lk*(hx[0]-hx[1]+ld*cos(phi)); <-- old formula, probably equivalent
    forcex[1]       =   -forcex[0];
    forcey[0]       =   lk * stretch * sin(phi); //-lk*(hy[0]-hy[1]+ld*sin(phi)); <-- old formula, probably equivalent
    forcey[1]       =   -forcey[0];

    force_par[0]    =   forcex[0]*actin_network->get_direction(aindex[0])[0] + forcey[0]*actin_network->get_direction(aindex[0])[1];
    force_perp[0]   =   -forcex[0]*actin_network->get_direction(aindex[0])[1] + forcey[0]*actin_network->get_direction(aindex[0])[0];
    force_par[1]    =   forcex[1]*actin_network->get_direction(aindex[1])[0] + forcey[1]*actin_network->get_direction(aindex[1])[1];
    force_perp[1]   =   -forcex[1]*actin_network->get_direction(aindex[1])[1] + forcey[1]*actin_network->get_direction(aindex[1])[0];

    torque[0]       =   cross(hx[0]-actin_network->get_position(aindex[0])[0],hy[0]-actin_network->get_position(aindex[0])[1],forcex[0],forcey[0]);
    torque[1]       =   cross(hx[1]-actin_network->get_position(aindex[1])[0],hy[1]-actin_network->get_position(aindex[1])[1],forcex[1],forcey[1]);

    actin_network->update_forces(aindex[0],force_par[0],force_perp[0],torque[0]);
    actin_network->update_forces(aindex[1],force_par[1],force_perp[1],torque[1]);
}
