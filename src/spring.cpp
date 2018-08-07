/*------------------------------------------------------------------
 spring.cpp : object describing a hookean spring that connects beads
            in a worm-like chain. 
 
 Copyright (C) 2016 
 Created by: Simon Freedman, Shiladitya Banerjee, Glen Hocky, Aaron Dinner
 Contact: dinner@uchicago.edu

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version. See ../LICENSE for details. 
-------------------------------------------------------------------*/

#include "spring.h"
#include "globals.h"
#include "filament.h"

spring::spring(){ }

spring::spring(double len, double stretching_stiffness, double max_ext_ratio, filament* f, 
        array<int, 2> myaindex, array<double, 2> myfov, array<int, 2> mynq)
{
    kl      = stretching_stiffness;
    l0      = len;
    fil     = f;
    aindex  = myaindex;
    fov     = myfov;
    nq      = mynq;
    half_nq = {{nq[0]/2, nq[1]/2}};

    max_ext = max_ext_ratio * l0;
    eps_ext = 0.01*max_ext;
    
    hx = {{0,0}};
    hy = {{0,0}};

    force = {{0,0}};
    intpoint = {{0,0}};
    point = {{0,0}}; 
    
    llen = l0;
    llensq = l0*l0;
}

spring::~spring(){ 
};

array<double,2> spring::get_hx(){
    return hx;
}

array<double,2> spring::get_hy(){
    return hy;
}

// stepping kinetics

void spring::step(string bc, double shear_dist)
{
    hx[0] = fil->get_bead(aindex[0])->get_xcm();
    hx[1] = fil->get_bead(aindex[1])->get_xcm();
    hy[0] = fil->get_bead(aindex[0])->get_ycm();
    hy[1] = fil->get_bead(aindex[1])->get_ycm();

    disp   = rij_bc(bc, hx[1]-hx[0], hy[1]-hy[0], fov[0], fov[1], shear_dist); 
    llensq = disp[0]*disp[0] + disp[1]*disp[1];
    llen   = sqrt(llensq);

    if (llen != 0)
        direc = {{disp[0]/llen, disp[1]/llen}};
    else
        direc = {{0, 0}};

}

void spring::update_force(string bc, double shear_dist)
{
    double kf = kl*(llen-l0);
    force = {{kf*direc[0], kf*direc[1]}};
}

/* Taken from hsieh, jain, larson, jcp 2006; eqn (5)
 * Adapted by placing a cutoff, similar to how it's done in LAMMPS src/bond_fene.cpp*/
void spring::update_force_fraenkel_fene(string bc, double shear_dist)
{
    double ext = abs(l0 - llen);
    double scaled_ext, klp;
    if (max_ext - ext > eps_ext ){
        scaled_ext = ext/max_ext;
    }
    else{
        scaled_ext = (max_ext - eps_ext)/max_ext;
    }

    klp = kl/(1-scaled_ext*scaled_ext)*(llen-l0);
    force = {{klp*direc[0], klp*direc[1]}};

}

void spring::update_force_marko_siggia(string bc, double shear_dist, double kToverLp)
{
    double xrat = llen/l0, yrat = llen/l0;
    if (xrat != xrat || xrat == 1) xrat = 0;
    if (yrat != yrat || yrat == 1) yrat = 0;
    force = {{kToverLp*(0.25/((1-xrat)*(1-xrat))-0.25+xrat), kToverLp*(0.25/((1-yrat)*(1-yrat))-0.25+yrat)}};  
}

array<double,2> spring::get_force()
{
    return force;
}

array<double,2> spring::get_disp()
{
    return disp;
}

array<double,2> spring::get_neg_disp()
{
    return {{-disp[0], -disp[1]}};
}

void spring::filament_update()
{
    fil->update_forces(aindex[0],  force[0],  force[1]);
    fil->update_forces(aindex[1], -force[0], -force[1]);

}

double spring::get_kl(){
    return kl;
}

void spring::set_l0(double myl0){
    l0 = myl0;
}

double spring::get_l0(){
    return l0;
}

double spring::get_fene_ext(){
    return max_ext/l0;
}

double spring::get_length(){
    return llen;
}

double spring::get_length_sq(){
    return llensq;
}

std::string spring::write(string bc, double shear_dist){
    return "\n" + std::to_string(hx[0]) + "\t" + std::to_string(hy[0]) + "\t" + std::to_string(disp[0]) + "\t" 
        + std::to_string(disp[1]);
}

std::string spring::to_string(){
    
    char buffer [100];
    sprintf(buffer, "aindex[0] = %d;\t aindex[1] = %d;\t kl = %f;\t l0 = %f\nfilament : \n",
                        aindex[0], aindex[1], kl, l0);

    return buffer + fil->to_string();

}

bool spring::operator==(const spring& that) 
{
    /*Note: you can't compare the filament objects because that will lead to infinite recursion;
     * this function requires the filament poiner to be identical to evaluate to true*/
    return (this->aindex[0] == that.aindex[0] && this->aindex[1] == that.aindex[1] &&
            this->kl == that.kl && 
            this->l0 == that.l0 && this->fil == that.fil);
}

bool spring::is_similar(const spring& that) 
{
    
    /* Same as ==; but doesn't compare the filament pointer*/

    return (this->aindex[0] == that.aindex[0] && this->aindex[1] == that.aindex[1] &&
            this->kl == that.kl &&
            this->l0 == that.l0);
}

// Updates all derived quantities of a monomer
void spring::quad_update(string bc, double delrx){
    
    //quadrant numbers crossed by the bead in x-direction
    quad.clear();
    int xlower, xupper, ylower, yupper;
    
    if (disp[0] >= 0)
    {
        xlower = coord2quad_floor(fov[0], nq[0], hx[0]);
        xupper = coord2quad_ceil(fov[0], nq[0], hx[1]);
    }
    else
    {
        xlower = coord2quad_floor(fov[0], nq[0], hx[1]);
        xupper = coord2quad_ceil(fov[0], nq[0], hx[0]);
    };
    
    if (disp[1] >=0 )
    {
        ylower = coord2quad_floor(fov[1], nq[1], hy[0]);
        yupper = coord2quad_ceil(fov[1], nq[1], hy[1]);
    }
    else
    {
        ylower = coord2quad_floor(fov[1], nq[1], hy[1]);
        yupper = coord2quad_ceil(fov[1], nq[1], hy[0]);
    };

    vector<int> xcoords = range_bc(bc, delrx, 0, nq[0], xlower, xupper);
    vector<int> ycoords = range_bc(bc, delrx, 0, nq[1], ylower, yupper);
    for(int xcoord : xcoords)
        for(int ycoord : ycoords){
            quad.push_back({{xcoord, ycoord}});
        }
}

//shortest(perpendicular) distance between an arbitrary point and the spring
//SO : 849211
double spring::get_distance_sq(string bc, double delrx, double xp, double yp)
{
    array<double, 2> dr = rij_bc(bc, intpoint[0]-xp, intpoint[1]-yp, fov[0], fov[1], delrx);
    return dr[0]*dr[0] + dr[1]*dr[1];
}

array<double,2> spring::get_intpoint()
{
    return intpoint;
}

void spring::calc_intpoint(string bc, double delrx, double xp, double yp)
{
    double l2 = disp[0]*disp[0]+disp[1]*disp[1];
    if (l2==0){
        intpoint = {{hx[0], hy[0]}};
    }else{
        //Consider the line extending the spring, parameterized as h0 + tp ( h1 - h0 )
        //tp = projection of (xp, yp) onto the line
        double tp=dot_bc(bc, xp-hx[0], yp-hy[0], hx[1]-hx[0], hy[1]-hy[0], fov[0], fov[1], delrx)/l2;
        if (tp<0){ 
            intpoint = {{hx[0], hy[0]}};
        }
        else if(tp>1.0){
            intpoint = {{hx[1], hy[1]}};
        }
        else{
            array<double, 2> proj   = {{hx[0] + tp*disp[0], hy[0] + tp*disp[1]}};
            intpoint                = pos_bc(bc, delrx, 0, fov, {{0,0}}, proj); //velocity and dt are 0 since not relevant
        }
    }
}

double spring::get_r_c(string bc, double delrx, double x, double y)
{
    double l2 = disp[0]*disp[0] + disp[1]*disp[1]; 
    double r_c; 
    double dx, dy; 
    array <double, 2> pos; 
    array <double, 2> proj; 
	
    if(l2 == 0)
    {
        point = {{hx[0], hy[0]}};
	pos = {{x,y}};
        dx = pos[0] - point[0];
        dy = pos[1] - point[1];
        r_c = dist_bc(bc, dx, dy, fov[0], fov[1], delrx);  
    }
    else
    {
        double tp = dot_bc(bc, x-hx[0], y-hy[0], hx[1]-hx[0], hy[1]-hy[0], fov[0], fov[1], delrx)/l2;

        if(tp < 0)
        {
            point = {{hx[0], hy[0]}};
            pos = {{x,y}};
            dx = pos[0] - point[0];
            dy = pos[1] - point[1];
            r_c = dist_bc(bc, dx, dy, fov[0], fov[1], delrx);
        }
        else if(tp > 1.0)
        {
            point = {{hx[1], hy[1]}};
            pos = {{x,y}};
            dx = pos[0] - point[0];
            dy = pos[1] - point[1];
            r_c = dist_bc(bc, dx, dy, fov[0], fov[1], delrx);
        }
        else
        {
            proj = {{hx[0] + tp*disp[0], hy[0] + tp*disp[1]}};
            point = pos_bc(bc, delrx, 0, fov, {{0,0}}, proj);
            pos = {{x,y}}; 
            dx = pos[0] - point[0];
            dy = pos[1] - point[1];
            r_c = dist_bc(bc, dx, dy, fov[0], fov[1], delrx);
        }

        //cout << "r_c: " << r_c << endl; 
	//cout << "tp: " << tp << endl; 
	//return r_c; 
    }
    return r_c; 
}

bool spring::get_line_intersect(string bc, double delrx, spring *l2)
{
    //Reference to Stack Overflow entry by iMalc on Feb 10, 2013
    //Web Address: https://stackoverflow.com/questions/563198/how-do-you-detect-where-two-line-segments-intersect

    //double dx1, dx2, dy1, dy2, dx12, 
    double dx12, dy12, denom, s_num, t_num;
    array <double,2> disp1, disp2, disp12, hx2, hy2;  
    bool denomPos; 
 
    //dx1 = hx[1]-hx[0];  
    //dy1 = hy[1]-hy[0]; 
    //dx2 = hx2[1]-hx2[0]; 
    //dy2 = hy2[1]-hy2[0]; 

    disp1 = this->get_disp(); 
    disp2 = l2->get_disp(); 

    //disp1 = rij_bc(bc, dx1, dy1, fov[0], fov[1], delrx); 
    //disp2 = rij_bc(bc, dx2, dy2, fov[0], fov[1], delrx); 

    hx2 = l2->get_hx();  
    hy2 = l2->get_hy(); 

    dx12 = hx[0]-hx2[0]; 
    dy12 = hy[0]-hy2[0]; 

    disp12 = rij_bc(bc, dx12, dy12, fov[0], fov[1], delrx); 

    denom  = (disp1[0]*disp2[1] - disp1[1]*disp2[0]); 
    if(denom == 0){return false;}
    denomPos = denom > 0;

    s_num = disp1[0]*disp12[1] - disp1[1]*disp12[0];   
    t_num = disp2[0]*disp12[1] - disp2[1]*disp12[0]; 

    if((s_num < 0) == denomPos){ return false; } 
    if((t_num < 0) == denomPos){ return false; }

    if(((s_num > denom) == denomPos) || ((t_num > denom) == denomPos)){ return false; }

    //Else Collision have been detected, the filaments do intersect!
    else{ return true; } 
     
}

array <double, 2> spring::get_point()
{
    //this->calc_r_c(bc,delrx,x,y);
    return point; 
}

vector<array<int, 2> > spring::get_quadrants()
{
    return quad;
}

array<double, 2> spring::get_direction()
{
    return direc;
}

double spring::get_stretching_energy(){
    return (force[0]*force[0]+force[1]*force[1])/(2*kl);
}

double spring::get_stretching_energy_fene(string bc, double shear_dist)
{
    double ext = abs(l0 - llen);
    
    if (max_ext - ext > eps_ext )
        return -0.5*kl*max_ext*max_ext*log(1-(ext/max_ext)*(ext/max_ext));
    else
        return 0.25*kl*ext*ext*(max_ext/eps_ext);
    
}

void spring::set_aindex(array<int,2> aidx)
{
    aindex = aidx;
}


double spring::get_max_ext()
{
    return max_ext;
}

array<int, 2> spring::get_aindex(){
    return aindex;
}

// functions for growing
void spring::inc_aindex()
{
    aindex = {{aindex[0]+1, aindex[1]+1}};
}

void spring::add_mot(motor * mot, int hd)
{
    mots[mot] = hd;
}

void spring::remove_mot(motor * mot)
{
    //cout<<"\nDEBUG: trying to remove motor";
///    cout<<"\nDEBUG: removing mot "<<mot<<" from mots on spring "<<this;
    mots.erase(mot);
}

map<motor *, int> & spring::get_mots()
{
    map<motor *, int> &ptr = mots;
    return ptr;
}
/*
int spring::get_n_mots(){
    return int(mots.size());
}

motor * spring::get_mot(int i){
    return mots[i];
}

int spring::get_mot_hd(int i){
    return mot_hds[i];
}*/
