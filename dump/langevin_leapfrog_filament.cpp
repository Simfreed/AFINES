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

langevin_leapfrog_filament::~langevin_leapfrog_filament(){
    cout<<"DELETING langevin_leapfrog_FILAMENT\n";
    for (unsigned int i = 0; i < actins.size(); i++) delete actins[i];
    for (unsigned int i = 0; i < links.size(); i++) delete links[i];
    actins.clear();
    links.clear();
}

langevin_leapfrog_filament::langevin_leapfrog_filament(vector<actin *> actinvec, array<double, 2> myfov, array<int, 2> mynq, double linkLength, 
        double stretching_stiffness, double bending_stiffness, 
        double deltat, double temp, double frac_force, double g, string bdcnd) : filament(actinvec, myfov, mynq, linkLength, stretching_stiffness, bending_stiffness, deltat, temp, frac_force, g, bdcnd)
{ 
    mass = actin_mass_density*linkLength;
    double fric_coef;
    if (actins.size() > 0) fric_coef = actins[0]->get_friction() / mass;
    else fric_coef = 0; 

    //Langevin Leap Frog coefficients, taken from Chodera et al
    a =          exp(-fric_coef*dt/2);
    b =     (1 - exp(-fric_coef*dt/2))/(fric_coef*dt);
    c = sqrt(1 - exp(-fric_coef*dt));
}

langevin_leapfrog_filament::langevin_leapfrog_filament(array<double, 3> startpos, int nactin, array<double, 2> myfov, array<int, 2> mynq, double visc, 
        double deltat, double temp, bool isStraight, double actinRadius, double linkLength, double stretching_stiffness,
        double bending_stiffness, double frac_force, string bdcnd):
    filament(startpos, nactin, myfov, mynq, visc, deltat, temp, isStraight, actinRadius, linkLength, stretching_stiffness, bending_stiffness,
            frac_force, bdcnd) 
{ 
    mass = actin_mass_density*linkLength;
    double fric_coef;
    if (actins.size() > 0) fric_coef = actins[0]->get_friction() / mass;
    else fric_coef = 0; 

    //Langevin Leap Frog coefficients, taken from Chodera et al
    a =          exp(-fric_coef*dt/2);
    b =     (1 - exp(-fric_coef*dt/2))/(fric_coef*dt);
    c = sqrt(1 - exp(-fric_coef*dt));
    cout<<"\nDEBUG: a = "<<a<<"\tb = "<<b<<"\tc = "<<c<<"\n";
}

langevin_leapfrog_filament::langevin_leapfrog_filament(array<double, 2> myfov, array<int, 2> mynq, double deltat, double temp, double shear, 
            double frac, double bending_stiffness, string bndcnd):
    filament(myfov, mynq, deltat, temp, shear, frac, bending_stiffness, bndcnd){ 
        a = b = c = mass = 0;
}

langevin_leapfrog_filament::langevin_leapfrog_filament():filament(){ a = b = c = mass = 0;}

bool langevin_leapfrog_filament::operator==(const langevin_leapfrog_filament& that){
    
    if (actins.size() != that.actins.size() || links.size() != that.links.size())
        return false;

    for (unsigned int i = 0; i < actins.size(); i++)
        if (!(*(actins[i]) == *(that.actins[i])))
            return false;
    
    for (unsigned int i = 0; i < links.size(); i++)
        if (!(links[i]->is_similar(*(that.links[i]))))
            return false;

    return (this->fov[0] == that.fov[0] && this->fov[1] == that.fov[1] && 
            this->nq[0] == that.nq[0] && this->nq[1] == that.nq[1] &&
            this->gamma == that.gamma && this->temperature == that.temperature &&
            this->dt == that.dt && this->fracture_force == that.fracture_force &&
            this->mass == that.mass);

}

void langevin_leapfrog_filament::set_mass(double m)
{
    mass = m;
}

void langevin_leapfrog_filament::reset_velocity()
{
    for (unsigned int i = 0; i < actins.size(); i++)
        actins[i]->reset_velocity();
}

void langevin_leapfrog_filament::update_velocity_brownian()
{
    for (unsigned int i = 0; i < actins.size(); i++){
        actins[i]->update_velocity( c*sqrt(temperature/mass)*rng_n(0,1), 
                                    c*sqrt(temperature/mass)*rng_n(0,1));
    }
}

void langevin_leapfrog_filament::update_velocity_drag()
{
    for (unsigned int i = 0; i < actins.size(); i++){
        actins[i]->update_velocity( a*actins[i]->get_velocity()[0],
                                    a*actins[i]->get_velocity()[1] );
    }
}

void langevin_leapfrog_filament::update_velocity_int_forces()
{
    for (unsigned int i = 0; i < actins.size(); i++){
        actins[i]->update_velocity( b*dt*actins[i]->get_force()[0]/mass,
                                    b*dt*actins[i]->get_force()[1]/mass);
    }
}
    
void langevin_leapfrog_filament::update_positions(double t)
{
    array<double, 2> newpos;
    kinetic_energy = 0;

    for (unsigned int i = 0; i < actins.size(); i++){

        kinetic_energy += actins[i]->get_vsquared();
        newpos = boundary_check(i, actins[i]->get_xcm() + actins[i]->get_velocity()[0]*dt, 
                                   actins[i]->get_ycm() + actins[i]->get_velocity()[1]*dt); 
        
        actins[i]->set_xcm(newpos[0]);
        actins[i]->set_ycm(newpos[1]);
        actins[i]->reset_force(); 

    }
    kinetic_energy *= mass*actins.size()/2;
    
    for (unsigned int i = 0; i < links.size(); i++)
        links[i]->step();

    
}

vector<langevin_leapfrog_filament *> langevin_leapfrog_filament::update_stretching(double t)
{
    vector<langevin_leapfrog_filament *> newfilaments;

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

vector<langevin_leapfrog_filament *> langevin_leapfrog_filament::fracture(int node){

    vector<langevin_leapfrog_filament *> newfilaments;
    cout<<"\n\tDEBUG: fracturing at node "<<node;
    
    if(links.size() == 0)
        return newfilaments;

    vector<actin *> lower_half = this->get_actins(0, node+1);
    vector<actin *> upper_half = this->get_actins(node+1, actins.size());

    if (lower_half.size() > 0)
        newfilaments.push_back(
                new langevin_leapfrog_filament(lower_half, fov, nq, links[0]->get_length(), links[0]->get_kl(), kb, 
                    dt, temperature, fracture_force, gamma, BC));
    if (upper_half.size() > 0)
        newfilaments.push_back(
                new langevin_leapfrog_filament(upper_half, fov, nq, links[0]->get_length(), links[0]->get_kl(), kb, 
                    dt, temperature, fracture_force, gamma, BC));

    for (int i = 0; i < (int)(lower_half.size()); i++) delete lower_half[i];
    for (int i = 0; i < (int)(upper_half.size()); i++) delete upper_half[i];
    
    lower_half.clear();
    upper_half.clear();
    
    return newfilaments;

}

