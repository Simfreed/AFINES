/*
 *  actin.cpp
 *  
 *
 *  Created by Shiladitya Banerjee on 9/3/13.
 *  Adapted into spherical setting by Simon Freedman, 3/2015
 *  Copyright 2013 University of Chicago. All rights reserved.
 *
 */

#include "actin.h"
#include "globals.h"

//actin filament class
actin::actin(){}

actin::actin(double xcm, double ycm, double len, double fovx, double fovy, int nx, int ny, double vis)
{
    //now i will make these spherical
    x=xcm;
    y=ycm;
    ld=len; //radius
    a_vis=vis;
    friction = 4*pi*a_vis*ld;
    
    //fov = new array<double, 2>();
    //nq = new array<int, 2>();
    fov[0] = fovx;
    fov[1] = fovy;
    nq[0]  = nx;
    nq[1]  = ny;
    
    forces[0] = 0;
    forces[1] = 0;
}

actin::actin(const actin& other){
    
    //cout<<"\nDEBUG: calling copy constructor"; 
    x = other.x;
    y = other.y;
    ld = other.ld;
    a_vis = other.a_vis;
    friction = other.friction;

    fov[0] = other.fov[0];
    fov[1] = other.fov[1];
    nq[0] = other.nq[0];
    nq[1] = other.nq[1];

    forces[0] = other.forces[0];
    forces[1] = other.forces[1];
    
}

actin::~actin(){ 
    //cout<<"DELETING ACTIN\n";
};


array<double,2> actin::get_forces()
{
    return forces;
}

array<double,2> actin::get_fov(){
    return fov;
}

array<int, 2> actin::get_nq(){
    return nq;
}


double actin::get_length()
{
    return ld;
}

void actin::update_force(double f1, double f2)
{
    forces[0]+=f1;
    forces[1]+=f2;
}

double actin::get_xcm()
{
    return x;
}

double actin::get_ycm()
{
    return y;
}

void actin::set_xcm(double xcm)
{
    x = xcm;
}

void actin::set_ycm(double ycm)
{
    y = ycm;
}

bool actin::operator==(const actin& that) 
{
    double err = eps; 
    return (close(x , that.x , err) && close(y , that.y , err) &&
            close(ld , that.ld , err) &&
            close(a_vis , that.a_vis , err) && close(forces[0] , that.forces[0] , err) &&
            close(forces[1] , that.forces[1] , err)
           );
}

string actin::write()
{
    return std::to_string(x) + "\t" + std::to_string(y);
}

string actin::to_string()
{
    return "x : " + std::to_string(x) + "\ty : " + std::to_string(y) +
           "\tld : " + std::to_string(ld) + "\ta_vis : "+ std::to_string(a_vis) + 
           "\tforces[0] : " + std::to_string(forces[0]) + "\tforces[1] : "+ std::to_string(forces[1]) + "\n";
 
}

double actin::get_viscosity(){
    return a_vis;
}
