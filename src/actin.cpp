/*------------------------------------------------------------------
 actin.cpp : object describing a circular bead
 
 Copyright (C) 2016 
 Created by: Simon Freedman, Shiladitya Banerjee, Glen Hocky, Aaron Dinner
 Contact: dinner@uchicago.edu

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version. See ../LICENSE for details. 
-------------------------------------------------------------------*/

#include "actin.h"
#include "globals.h"

actin::actin(){}

actin::actin(double xcm, double ycm, double len, double vis)
{
    //now i will make these spherical
    x=xcm;
    y=ycm;
    ld=len; //radius
    a_vis=vis;
    friction = 6*pi*a_vis*ld;
   
    force = {0,0};
    velocity = {0,0};
}

actin::actin(const actin& other){
    
    //cout<<"\nDEBUG: calling copy constructor"; 
    x = other.x;
    y = other.y;
    ld = other.ld;
    a_vis = other.a_vis;
    friction = other.friction;
    force = other.force;
    velocity = other.velocity;
}

actin::~actin(){ 
    //cout<<"DELETING ACTIN\n";
};


array<double,2> actin::get_force()
{
    return force;
}

array<double,2> actin::get_velocity()
{
    return velocity;
}

double actin::get_vsquared()
{
    return velocity[0]*velocity[0] + velocity[1]*velocity[1];
}

double actin::get_length()
{
    return ld;
}

void actin::update_velocity(double v1, double v2)
{
    velocity[0]+=v1;
    velocity[1]+=v2;
}

void actin::update_force(double f1, double f2)
{
    if(f1 == f1 && f2 == f2 && std::isfinite(f1) && std::isfinite(f2)){
        force[0]+=f1;
        force[1]+=f2;
    }else{
        cout<<"\nENCOUNTERED INFINITE FORCE; PROGRAM ABORTING\n";
        abort();
    }
}

void actin::reset_velocity()
{
    velocity[0] = 0;
    velocity[1] = 0;
}

void actin::reset_force()
{
    force[0] = 0;
    force[1] = 0;
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
    return (close( this->x , that.x , err) && close( this->y , that.y , err) &&
            close( this->ld , that.ld , err) &&
            close( this->a_vis , that.a_vis , err) && close( this->force[0] , that.force[0] , err) &&
            close( this->force[1] , that.force[1] , err)
           );
}

string actin::write()
{
    return "\n" + std::to_string(x) + "\t" + std::to_string(y) + "\t" + std::to_string(ld);
}

string actin::to_string()
{
    return "x : " + std::to_string(x) + "\ty : " + std::to_string(y) +
           "\tld : " + std::to_string(ld) + "\ta_vis : "+ std::to_string(a_vis) + 
           "\tforce[0] : " + std::to_string(force[0]) + "\tforce[1] : "+ std::to_string(force[1]) + "\n";
 
}

double actin::get_viscosity(){
    return a_vis;
}

double actin::get_friction(){
    return friction;
}

double actin::get_ld(){
    return ld;
}
