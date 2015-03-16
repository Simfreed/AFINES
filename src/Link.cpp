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

Link::Link(double len, double stretching_stiffness, filament* f, 
        array<int, 2> myaindex, array<double, 2> myfov, array<int, 2> mynq)
{
    kl      = stretching_stiffness;
    l0      = len;
    fil     = f;
    aindex  = myaindex;
    fov     = myfov;
    nq      = mynq;

    hx = {0,0};
    hy = {0,0};

    this->step();
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

void Link::set_aindex1(int i){
    aindex[1] = i;
    this->step();
}

// stepping kinetics
void Link::step()
{
    hx[0] = fil->get_actin(aindex[0])->get_xcm();
    hx[1] = fil->get_actin(aindex[1])->get_xcm();
    hy[0] = fil->get_actin(aindex[0])->get_ycm();
    hy[1] = fil->get_actin(aindex[1])->get_ycm();

    xcm = (hx[0]+hx[1])/2.0;
    ycm = (hy[0]+hy[1])/2.0;
    phi=atan2(hy[1]-hy[0],hx[1]-hx[0]);

}

double Link::get_stretch_force(){
    //cout<<"\nDEBUG:Stretch force = "<<kl*(dis_points(hx[0],hy[0],hx[1],hy[1])-l0);
    if (kl!=kl) cout<<"\nDEBUG: kl is inf";
    return kl * (dis_points(hx[0],hy[0],hx[1],hy[1])-l0);
}

void Link::filament_update()
{

    double force_stretch = this->get_stretch_force();
    double fx0, fy0, fx1, fy1;
    
    fx0 =  -force_stretch * cos(phi); 
    fy0 =  -force_stretch * sin(phi); 
    fx1 =  -fx0;
    fy1 =  -fy0;
    fil->update_forces(aindex[0], fx0, fy0);
    fil->update_forces(aindex[1], fx1, fy1);

}

double Link::get_kl(){
    return kl;
}

double Link::get_l0(){
    return l0;
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
    return dis_points(hx[0], hy[0], hx[1], hy[1]);
}

std::string Link::write(){
    return std::to_string(hx[0]) + "\t" + std::to_string(hy[0]) + "\t" + std::to_string(hx[1]-hx[0]) + "\t" 
        + std::to_string(hy[1]-hy[0]) + "\n";
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
void Link::quad_update(){
    
    //quadrant numbers crossed by the actin in x-direction
    quad.clear();
    int xlower, xupper, ylower, yupper;
    
    if(hx[0] <= hx[1])
    {
        xlower = int(floor(hx[0]/fov[0]*nq[0]));
        xupper = int(ceil( hx[1]/fov[0]*nq[0]));
    }
    else
    {
        xlower = int(floor(hx[1]/fov[0]*nq[0]));
        xupper = int(ceil( hx[0]/fov[0]*nq[0]));
    };
    
    if(hy[0] <= hy[1])
    {
        ylower = int(floor(hy[0]/fov[1]*nq[1]));
        yupper = int(ceil( hy[1]/fov[1]*nq[1]));
    }
    else
    {
        ylower = int(floor(hy[1]/fov[1]*nq[1]));
        yupper = int(ceil( hy[0]/fov[1]*nq[1]));
    };
    
    for(int xcoord = xlower; xcoord < xupper; xcoord++)
        for(int ycoord = ylower; ycoord < yupper; ycoord++)
            quad.push_back({xcoord, ycoord});
        
}

//shortest(perpendicular) distance between an arbitray point and the Link
double Link::get_distance(double xp, double yp)
{
    double l2=pow(dis_points(hx[0],hy[0],hx[1],hy[1]),2);
    if (l2==0) {
        return dis_points(xp,yp,hx[0],hy[0]);
    }
    double tp=dot(xp-hx[0],yp-hy[0],hx[1]-hx[0],hy[1]-hy[0])/l2;
    if (tp<0) {
        return dis_points(xp,yp,hx[0],hy[0]);
    }
    else if(tp>1.0){
        return dis_points(xp,yp,hx[1],hy[1]);
    }
    else{
        double px=hx[0]+tp*(hx[1]-hx[0]);
        double py=hy[0]+tp*(hy[1]-hy[0]);
        return dis_points(xp,yp,px,py);
    }
}

//closest point on the link to point (xp, yp)
array<double,2> Link::get_intpoint(double xp, double yp)
{
    array<double,2> coordinates;
    double l2 = pow(dis_points(hx[0], hy[0], hx[1], hy[1]) , 2);
    if (l2==0) {
        coordinates[0]=hx[0];
        coordinates[1]=hy[1];
    }
    double tp=dot(xp-hx[0],yp-hy[0],hx[1]-hx[0],hy[1]-hy[0])/l2;
    if (tp<0) {
        coordinates[0]=hx[0];
        coordinates[1]=hy[0];
    }
    else if(tp>1.0){
        coordinates[0]=hx[1];
        coordinates[1]=hy[1];
    }
    else{
        coordinates[0]=hx[0]+tp*(hx[1]-hx[0]);
        coordinates[1]=hy[0]+tp*(hy[1]-hy[0]);
    }
    return coordinates;
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
