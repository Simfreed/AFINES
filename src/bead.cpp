/*------------------------------------------------------------------
 bead.cpp : object describing a circular bead
 
 Copyright (C) 2016 
 Created by: Simon Freedman, Shiladitya Banerjee, Glen Hocky, Aaron Dinner
 Contact: dinner@uchicago.edu

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version. See ../LICENSE for details. 
-------------------------------------------------------------------*/

#include "bead.h"
#include "globals.h"

bead::bead(){}

bead::bead(double xcm, double ycm, double len, double vis)
{
    x=xcm;
    y=ycm;
    rad=len; //radius
    visc=vis;
    friction = 6*pi*visc*rad;
   
    force = {{0,0}};
}

bead::bead(const bead& other){
    
    x = other.x;
    y = other.y;
    rad = other.rad;
    visc = other.visc;
    friction = other.friction;
    force = other.force;
}

bead::~bead(){ 
};


array<double,2> bead::get_force()
{
    return force;
}

double bead::get_length()
{
    return rad;
}

void bead::update_force(double f1, double f2)
{
    if(f1 == f1 && f2 == f2 && std::isfinite(f1) && std::isfinite(f2)){
        force[0]+=f1;
        force[1]+=f2;
    }else{
        cout<<"\nENCOUNTERED INFINITE FORCE; PROGRAM ABORTING\n";
        abort();
    }
}

void bead::reset_force()
{
    force[0] = 0;
    force[1] = 0;
}

double bead::get_xcm()
{
    return x;
}

double bead::get_ycm()
{
    return y;
}

void bead::set_xcm(double xcm)
{
    x = xcm;
}

void bead::set_ycm(double ycm)
{
    y = ycm;
}

bool bead::operator==(const bead& that) 
{
    double err = eps; 
    return (close( this->x , that.x , err) && close( this->y , that.y , err) &&
            close( this->rad , that.rad , err) &&
            close( this->visc , that.visc , err) && close( this->force[0] , that.force[0] , err) &&
            close( this->force[1] , that.force[1] , err)
           );
}

string bead::write()
{
    return "\n" + std::to_string(x) + "\t" + std::to_string(y) + "\t" + std::to_string(rad);
}

string bead::to_string()
{
    return "x : " + std::to_string(x) + "\ty : " + std::to_string(y) +
           "\trad : " + std::to_string(rad) + "\tvisc : "+ std::to_string(visc) + 
           "\tforce[0] : " + std::to_string(force[0]) + "\tforce[1] : "+ std::to_string(force[1]) + "\n";
 
}

double bead::get_viscosity(){
    return visc;
}

double bead::get_friction(){
    return friction;
}
