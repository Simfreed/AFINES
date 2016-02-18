/*
 *  Link.cpp
 *  
 *
 *  Created by Simon Freedman on 9/10/2014
 *  Copyright 2014 University of Chicago. All rights reserved.
 *
 */

#include "Link.h"
#include "globals.h"
#include "filament.h"

Link::Link(){ }

Link::Link(double len, double stretching_stiffness, double max_ext_ratio, filament* f, 
        array<int, 2> myaindex, array<double, 2> myfov, array<int, 2> mynq)
{
    kl      = stretching_stiffness;
    l0      = len;
    fil     = f;
    aindex  = myaindex;
    fov     = myfov;
    nq      = mynq;

    max_ext = max_ext_ratio * l0;
    eps_ext = 0.01*max_ext;
    
    hx = {0,0};
    hy = {0,0};

    force = {0,0};
    //this->step();
}
Link::~Link(){ 
    //std::cout<<"DELETING LINK\n";
};

array<double,2> Link::get_hx(){
    return hx;
}

array<double,2> Link::get_hy(){
    return hy;
}

// stepping kinetics

void Link::step(string bc, double shear_dist)
{
    hx[0] = fil->get_actin(aindex[0])->get_xcm();
    hx[1] = fil->get_actin(aindex[1])->get_xcm();
    hy[0] = fil->get_actin(aindex[0])->get_ycm();
    hy[1] = fil->get_actin(aindex[1])->get_ycm();

    array<double, 2> cm = cm_bc(bc, {hx[0], hx[1]}, {hy[0], hy[1]}, fov[0], fov[1], shear_dist);
    xcm = cm[0];
    ycm = cm[1];
    
    disp = rij_bc(bc, hx[1]-hx[0], hy[1]-hy[0], fov[0], fov[1], shear_dist); 
    phi=atan2(disp[1],disp[0]);

}

void Link::update_force(string bc, double shear_dist)
{
    force = {kl*(disp[0]-l0*cos(phi)), kl*(disp[1]-l0*sin(phi))};
}

/* Taken from hsieh, jain, larson, jcp 2006; eqn (5)
 * Adapted by placing a cutoff, similar to how it's done in LAMMPS src/bond_fene.cpp*/
void Link::update_force_fraenkel_fene(string bc, double shear_dist)
{
    double ext = abs(l0 - hypot(disp[0], disp[1]));
    double scaled_ext, klp;
    if (max_ext - ext > eps_ext ){
        //cout<<"\nDEBUG: approaching cutoff";
        scaled_ext = ext/max_ext;
    }
    else{
        //cout<<"\nDEBUG: at cutoff";
        scaled_ext = (max_ext - eps_ext)/max_ext;
    }
    klp = kl/(1-scaled_ext*scaled_ext);
    force = {klp*(disp[0]-l0*cos(phi)), klp*(disp[1]-l0*sin(phi))};

}

void Link::update_force_marko_siggia(string bc, double shear_dist, double kToverLp)
{
    double xrat = disp[0]/(l0*cos(phi)), yrat = disp[1]/(l0*sin(phi));
    if (xrat != xrat || xrat == 1) xrat = 0;
    if (yrat != yrat || yrat == 1) yrat = 0;
    force = {kToverLp*(0.25/((1-xrat)*(1-xrat))-0.25+xrat), kToverLp*(0.25/((1-yrat)*(1-yrat))-0.25+yrat)};  
}

array<double,2> Link::get_force()
{
    return force;
}

array<double,2> Link::get_disp()
{
    return disp;
}

array<double,2> Link::get_neg_disp()
{
    return {-disp[0], -disp[1]};
}

void Link::filament_update()
{
    fil->update_forces(aindex[0],  force[0],  force[1]);
    fil->update_forces(aindex[1], -force[0], -force[1]);

}

double Link::get_kl(){
    return kl;
}

double Link::get_l0(){
    return l0;
}

double Link::get_fene_ext(){
    return max_ext/l0;
}

double Link::get_xcm(){
    return xcm;
}

double Link::get_ycm(){
    return ycm;
}

double Link::get_angle(){
    return phi;
}

double Link::get_length(){
    return l0; 
}

std::string Link::write(string bc, double shear_dist){
    return "\n" + std::to_string(hx[0]) + "\t" + std::to_string(hy[0]) + "\t" + std::to_string(disp[0]) + "\t" 
        + std::to_string(disp[1]);
}

std::string Link::to_string(){
    
    char buffer [100];
    sprintf(buffer, "aindex[0] = %d;\t aindex[1] = %d;\t kl = %f;\t l0 = %f\nfilament : \n",
                        aindex[0], aindex[1], kl, l0);

    return buffer + fil->to_string();

}

bool Link::operator==(const Link& that) 
{
    /*Note: you can't compare the filament objects because that will lead to infinite recursion;
     * this function requires the filament poiner to be identical to evaluate to true*/
    return (this->aindex[0] == that.aindex[0] && this->aindex[1] == that.aindex[1] &&
            this->kl == that.kl && 
            this->l0 == that.l0 && this->fil == that.fil);
}

bool Link::is_similar(const Link& that) 
{
    
    /* Same as ==; but doesn't compare the filament pointer*/

    return (this->aindex[0] == that.aindex[0] && this->aindex[1] == that.aindex[1] &&
            this->kl == that.kl &&
            this->l0 == that.l0);
}

//All these things swtich from being applicable to actin to being applicable to links:
// Updates all derived quantities of a monomer
void Link::quad_update(string bc, double delrx){
    
    //quadrant numbers crossed by the actin in x-direction
    quad.clear();
    int xlower, xupper, ylower, yupper;
    
    if(fabs(phi) < pi/2)//if(hx[0] <= hx[1])
    {
        xlower = int(round( hx[0]/fov[0]*nq[0]));
        xupper = int(round( hx[1]/fov[0]*nq[0]));
    }
    else
    {
        xlower = int(round( hx[1]/fov[0]*nq[0]));
        xupper = int(round( hx[0]/fov[0]*nq[0]));
    };
    
    if(phi >= 0) //hy[0] <= hy[1])
    {
        ylower = int(round( hy[0]/fov[1]*nq[1]));
        yupper = int(round( hy[1]/fov[1]*nq[1]));
    }
    else
    {
        ylower = int(round( hy[1]/fov[1]*nq[1]));
        yupper = int(round( hy[0]/fov[1]*nq[1]));
    };

    if (xlower == xupper) xupper++;
    if (ylower == yupper) yupper++;

    vector<int> xcoords = range_bc(bc, delrx, nq[0]/2, xlower, xupper);
    vector<int> ycoords = range_bc(bc, delrx, nq[1]/2, ylower, yupper);
    
    for(int xcoord : xcoords)
        for(int ycoord : ycoords){
            quad.push_back({xcoord, ycoord});
            //cout<<"\nDEBUG: quadrant : ("<<xcoord<<" , "<<ycoord<<")";
        }


}

//shortest(perpendicular) distance between an arbitrary point and the Link
//SO : 849211
double Link::get_distance(string bc, double delrx, double xp, double yp)
{
    array<double, 2> ip = this->get_intpoint(bc, delrx, xp, yp);
    return dist_bc(bc, ip[0]-xp, ip[1]-yp, fov[0], fov[1], delrx);
}

array<double,2> Link::get_intpoint(string bc, double delrx, double xp, double yp)
{
    array<double, 2> int_point; 
    double l2 = pow(dist_bc(bc, hx[1]-hx[0], hy[1]-hy[0], fov[0], fov[1], delrx),2);
    if (l2==0){
        int_point = {hx[0], hy[0]};
    }else{
        //Consider the line extending the link, parameterized as h0 + tp ( h1 - h0 )
        //tp = projection of (xp, yp) onto the line
        double tp=dot_bc(bc, xp-hx[0], yp-hy[0], hx[1]-hx[0], hy[1]-hy[0], fov[0], fov[1], delrx)/l2;

        if (tp<0){ 
            int_point = {hx[0], hy[0]};
        }
        else if(tp>1.0){
            int_point = {hx[1], hy[1]};
        }
        else{
            array<double, 2> proj   = {hx[0] + tp*disp[0], hy[0] + tp*disp[1]};
            int_point               = pos_bc(bc, delrx, 0, fov, {0,0}, proj); //velocity and dt are 0 since not relevant
            //cout<<"\nDEBUG: tp = "<<tp<<"; h0 = ("<<hx[0]<<", "<<hy[0]<<")\nproj = ("<<proj[0]<<", "<<proj[1]<<"); closest = ("<<closest[0]<<", "<<closest[1]<<")";
        }
    }
    return int_point;
}

double Link::get_int_angle(double xp, double yp)
{
    double angle;
    double xcor,ycor;
    double slope=(hy[1]-hy[0])/(hx[1]-hx[0]);
    double yintercept = ycm - slope * xcm;
    xcor=(slope*yp + xp - slope*yintercept)/(slope*slope + 1);
    ycor=(slope*slope*yp + slope*xp + yintercept)/(1 + slope*slope);
    angle=atan2((ycor-yp),(xcor-xp));
    return angle;
}

vector<array<int, 2> > Link::get_quadrants()
{
    return quad;
}

array<double, 2> Link::get_direction()
{
    return {cos(phi), sin(phi)};
}

double Link::get_stretching_energy(){
    return (force[0]*force[0]+force[1]*force[1])/(2*kl);
}

double Link::get_stretching_energy_fene(string bc, double shear_dist)
{
    double ext = abs(l0 - dist_bc(bc, hx[1]-hx[0], hy[1]-hy[0], fov[0], fov[1], shear_dist));
    
    if (max_ext - ext > eps_ext )
        return -0.5*kl*max_ext*max_ext*log(1-(ext/max_ext)*(ext/max_ext));
    else
        return 0.25*kl*ext*(max_ext/eps_ext);
    
}
