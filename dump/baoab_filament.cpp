/*
 *  filament.cpp
 *  
 *
 *  Created by Simon Freedman on 12/22/2014
 *  Copyright 2014 University of Chicago. All rights reserved.
 *
 */

#include "filament.h"
#include "actin.h"
#include "globals.h"

baoab_filament::~baoab_filament(){
    for (unsigned int i = 0; i < actins.size(); i++) delete actins[i];
    for (unsigned int i = 0; i < links.size(); i++) delete links[i];
    actins.clear();
    links.clear();
}

baoab_filament::baoab_filament(vector<actin *> actinvec, array<double, 2> myfov, array<int, 2> mynq, double linkLength, 
        double stretching_stiffness, double bending_stiffness, 
        double deltat, double temp, double frac_force, double g, string bdcnd) : filament(actinvec, myfov, mynq, linkLength, stretching_stiffness, bending_stiffness, deltat, temp, frac_force, g, bdcnd)
{ 
    mass = actin_mass_density*linkLength;
    double fric_coef = 0;
    if (actins.size() > 0) fric_coef = actins[0]->get_friction() / mass;
    a = exp(-fric_coef * dt);
    b = sqrt(1-a*a);
}
baoab_filament::baoab_filament(array<double, 3> startpos, int nactin, array<double, 2> myfov, array<int, 2> mynq, double visc, 
        double deltat, double temp, bool isStraight, double actinRadius, double linkLength, double stretching_stiffness,
        double bending_stiffness, double frac_force, string bdcnd):
    filament(startpos, nactin, myfov, mynq, visc, deltat, temp, isStraight, actinRadius, linkLength, stretching_stiffness, bending_stiffness,
            frac_force, bdcnd)  
{ 
    mass = actin_mass_density*linkLength;
    double fric_coef = 0;
    if (actins.size() > 0) fric_coef = actins[0]->get_friction() / mass;
    a = exp(-fric_coef * dt);
    b = sqrt(1-a*a);
}

/* BAOAB
 * vx = vx + Fx * dt /2
 * vy = vy + Fy * dt /2
 *
 * x = x + vx * dt /2
 * y = y + vy * dt /2
 *
 * vx = sqrt( temperature ) * rng_n(0,1)
 * vy = sqrt( temperature ) * rng_n(0,1)
 *
 * x = x + vx * dt /2
 * y = y + vy * dt /2
 *
 * vx = vx + Fx * dt /2
 * vy = vy + Fy * dt /2
 */

void baoab_filament::update_velocities_B()
{
    kinetic_energy = 0;  
    for (unsigned int i = 0; i < actins.size(); i++){
        actins[i]->update_velocity((actins[i]->get_force()[0])*dt/(2*mass), 
                                   (actins[i]->get_force()[1])*dt/(2*mass));
        
        kinetic_energy += actins[i]->get_vsquared();
    
    }
    kinetic_energy *= mass*actins.size()/2;
}

void baoab_filament::update_positions(double t)
{
    array<double, 2> newpos;
    dt = dt/2;
    
    for (unsigned int i = 0; i < actins.size(); i++){
        newpos = boundary_check(i, actins[i]->get_xcm() + actins[i]->get_velocity()[0]*dt, 
                                   actins[i]->get_ycm() + actins[i]->get_velocity()[1]*dt); 
        actins[i]->set_xcm(newpos[0]);
        actins[i]->set_ycm(newpos[1]);
        actins[i]->reset_force(); 
    }
    dt = 2*dt;
    
    for (unsigned int i = 0; i < links.size(); i++)
        links[i]->step();

    
}

void baoab_filament::update_velocities_O(double t)
{
    double T=temperature;
    array<double, 2> v;
    //if (t < dt*100000) T = 0;
    
    for (unsigned int i = 0; i < actins.size(); i++){
        v = actins[i]->get_velocity();
        actins[i]->reset_velocity();
        actins[i]->update_velocity(a * v[0] + b * sqrt(T/mass) * rng_n(0, 1), 
                                   a * v[1] + b * sqrt(T/mass) * rng_n(0, 1));
    }
}

vector<baoab_filament *> baoab_filament::update_stretching(double t)
{
    vector<baoab_filament *> newfilaments;
    
    if(links.size() == 0)
        return newfilaments;
   
    for (unsigned int i=0; i < links.size(); i++) {
        links[i]->update_force(BC, delrx);
        if (hypot(links[i]->get_force()[0], links[i]->get_force()[1]) > fracture_force){
            newfilaments = this->fracture(i);
            break;
        }
        else 
            links[i]->filament_update();
    }
    
    return newfilaments;
}

vector<baoab_filament *> baoab_filament::fracture(int node){

    vector<baoab_filament *> newfilaments;
    cout<<"\n\tDEBUG: fracturing at node "<<node;
    
    if(links.size() == 0)
        return newfilaments;

    vector<actin *> lower_half = this->get_actins(0, node+1);
    vector<actin *> upper_half = this->get_actins(node+1, actins.size());

    if (lower_half.size() > 0)
        newfilaments.push_back(
                new baoab_filament(lower_half, fov, nq, links[0]->get_length(), links[0]->get_kl(), kb, 
                    dt, temperature, fracture_force, gamma, BC));
    if (upper_half.size() > 0)
        newfilaments.push_back(
                new baoab_filament(upper_half, fov, nq, links[0]->get_length(), links[0]->get_kl(), kb, 
                    dt, temperature, fracture_force, gamma, BC));

    for (int i = 0; i < (int)(lower_half.size()); i++) delete lower_half[i];
    for (int i = 0; i < (int)(upper_half.size()); i++) delete upper_half[i];
    
    lower_half.clear();
    upper_half.clear();
    
    return newfilaments;

}

