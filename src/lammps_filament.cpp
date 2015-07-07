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

lammps_filament::~lammps_filament(){
    cout<<"DELETING LAMMPS_FILAMENT\n";
    for (unsigned int i = 0; i < actins.size(); i++) delete actins[i];
    for (unsigned int i = 0; i < links.size(); i++) delete links[i];
    actins.clear();
    links.clear();
}

lammps_filament::lammps_filament(vector<actin *> actinvec, array<double, 2> myfov, array<int, 2> mynq, double linkLength, 
        double stretching_stiffness, double bending_stiffness, 
        double deltat, double temp, double frac_force, double g, string bdcnd) : filament(actinvec, myfov, mynq, linkLength, stretching_stiffness, bending_stiffness, deltat, temp, frac_force, g, bdcnd)
{ 
    mass = actin_mass_density*linkLength;
}

lammps_filament::lammps_filament(array<double, 3> startpos, int nactin, array<double, 2> myfov, array<int, 2> mynq, double visc, 
        double deltat, double temp, bool isStraight, double actinRadius, double linkLength, double stretching_stiffness,
        double bending_stiffness, double frac_force, string bdcnd):
    filament(startpos, nactin, myfov, mynq, visc, deltat, temp, isStraight, actinRadius, linkLength, stretching_stiffness, bending_stiffness,
            frac_force, bdcnd) 
{ 
    mass = actin_mass_density*linkLength;
}

lammps_filament::lammps_filament(array<double, 2> myfov, array<int, 2> mynq, double deltat, double temp, double shear, 
            double frac, double bending_stiffness, string bndcnd):
    filament(myfov, mynq, deltat, temp, shear, frac, bending_stiffness, bndcnd){ mass = 0;}

lammps_filament::lammps_filament():filament(){ mass = 0;}

bool lammps_filament::operator==(const lammps_filament& that){
    
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

void lammps_filament::set_mass(double m)
{
    mass = m;
}

void lammps_filament::update_brownian()
{
    double gamma;
    for (unsigned int i = 0; i < actins.size(); i++){
        gamma = actins[i]->get_friction();  
        actins[i]->update_force(sqrt(2*temperature*gamma/dt)*rng_n(0,1), 
                                sqrt(2*temperature*gamma/dt)*rng_n(0,1));
    }
}

void lammps_filament::update_drag()
{
    double gamma;
    for (unsigned int i = 0; i < actins.size(); i++){
        gamma = -1*actins[i]->get_friction();  
        actins[i]->update_force( gamma*actins[i]->get_velocity()[0],
                                 gamma*actins[i]->get_velocity()[1] );
    }
}

void lammps_filament::update_positions(double t)
{
    double dtfm = dt/(2*mass);
    array<double, 2> newpos;
    kinetic_energy = 0;

    for (unsigned int i = 0; i < actins.size(); i++){

        actins[i]->update_velocity(actins[i]->get_force()[0]*dtfm,
                                   actins[i]->get_force()[1]*dtfm ); 
                    
        kinetic_energy += actins[i]->get_vsquared();
        
        newpos = boundary_check(i, t, actins[i]->get_velocity()[0], actins[i]->get_velocity()[1]); 
        
        actins[i]->set_xcm(newpos[0]);
        actins[i]->set_ycm(newpos[1]);
        actins[i]->reset_force(); 
        //actins[i]->reset_velocity(); 

    }
    
    for (unsigned int i = 0; i < links.size(); i++)
        links[i]->step();

    
}

vector<lammps_filament *> lammps_filament::update_stretching()
{
    vector<lammps_filament *> newfilaments;
    
    if(links.size() == 0)
        return newfilaments;
    
    for (unsigned int i=0; i < links.size(); i++) {
        if (fabs(links[i]->get_stretch_force()) > fracture_force){
            newfilaments = this->fracture(i);
            break;
        }
        else 
            links[i]->filament_update();
    }
    
    return newfilaments;
}

vector<lammps_filament *> lammps_filament::fracture(int node){

    vector<lammps_filament *> newfilaments;
    cout<<"\n\tDEBUG: fracturing at node "<<node;
    
    if(links.size() == 0)
        return newfilaments;

    vector<actin *> lower_half = this->get_actins(0, node+1);
    vector<actin *> upper_half = this->get_actins(node+1, actins.size());

    if (lower_half.size() > 0)
        newfilaments.push_back(
                new lammps_filament(lower_half, fov, nq, links[0]->get_length(), links[0]->get_kl(), kb, 
                    dt, temperature, fracture_force, gamma, BC));
    if (upper_half.size() > 0)
        newfilaments.push_back(
                new lammps_filament(upper_half, fov, nq, links[0]->get_length(), links[0]->get_kl(), kb, 
                    dt, temperature, fracture_force, gamma, BC));

    for (int i = 0; i < (int)(lower_half.size()); i++) delete lower_half[i];
    for (int i = 0; i < (int)(upper_half.size()); i++) delete upper_half[i];
    
    lower_half.clear();
    upper_half.clear();
    
    return newfilaments;

}

