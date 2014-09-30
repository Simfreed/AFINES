/*
 *  actin.cpp
 *  
 *
 *  Created by Shiladitya Banerjee on 9/3/13.
 *  Copyright 2013 University of Chicago. All rights reserved.
 *
 */

#include "actin.h"
#include "globals.h"

//actin filament class
actin::actin(double xcm, double ycm, double angle, double len, double fovx, double fovy, int nx, int ny, double vis)
{
    x=xcm;
    y=ycm;
    phi=angle;
    ld=len;
    diameter=ld/40;
    start[0]=x-ld*0.5*cos(phi);
    start[1]=y-ld*0.5*sin(phi);
    end[0]=x+ld*0.5*cos(phi);
    end[1]=y+ld*0.5*sin(phi);
    a_vis=vis;
    //unit vector
    e[0]=cos(phi);
    e[1]=sin(phi);
    //unit normal
    n[0] = -e[1];
    n[1] = e[0];
    //motor-induced forces
    forces[0]=0; //along the filament
    forces[1]=0; //perpendicular to the filament
    forces[2]=0; //torque

    //quadrant numbers crossed by the actin in x-direction
    quad.clear();
    tmp.clear();
    int lower_limit, upper_limit, index;
    if(start[0] <= end[0])
    {
        lower_limit = int(floor(start[0]/fovx*nx));
        if(lower_limit > 0){lower_limit--;};
        upper_limit = int(ceil(end[0]/fovx*nx));
        if(upper_limit < nx-1){upper_limit++;};

        for(index = lower_limit; index < upper_limit;index++){tmp.push_back(index);};
    }
    else
    {
        lower_limit = int(floor(end[0]/fovx*nx));
        if(lower_limit > 0){lower_limit--;};
        upper_limit = int(ceil(start[0]/fovx*nx));
        if(upper_limit < nx-1){upper_limit++;};
        for(index = lower_limit; index < upper_limit;index++){tmp.push_back(index);};
    };
    quad.push_back(tmp);

    //quadrant numbers crossed by the actin in y-direction
    tmp.clear();
    if(start[1] <= end[1])
    {
        lower_limit = int(floor(start[1]/fovy*ny));
        if(lower_limit > 0){lower_limit--;};
        upper_limit = int(ceil(end[1]/fovy*ny));
        if(upper_limit < ny-1){upper_limit++;};
        for(index = lower_limit; index < upper_limit;index++){tmp.push_back(index);}
    }
    else
    {
        lower_limit = int(floor(end[1]/fovy*ny));
        if(lower_limit > 0){lower_limit--;};
        upper_limit = int(ceil(start[1]/fovy*ny));
        if(upper_limit < ny-1){upper_limit++;};
        for(index = lower_limit; index < upper_limit;index++){tmp.push_back(index);}
    };
    quad.push_back(tmp);
}

actin::~actin(){ 
//    std::cout<<"DELETING ACTIN\n";
};

//shortest(perpendicular) distance between an arbitray point and the filament
double actin::get_distance(double xp, double yp)
{
    double l2=pow(dis_points(start[0],start[1],end[0],end[1]),2);
    if (l2==0) {
        return dis_points(xp,yp,start[0],start[1]);
    }
    double tp=dot(xp-start[0],yp-start[1],end[0]-start[0],end[1]-start[1])/l2;
    if (tp<0) {
        return dis_points(xp,yp,start[0],start[1]);
    }
    else if(tp>1.0){
        return dis_points(xp,yp,end[0],end[1]);
    }
    else{
        double px=start[0]+tp*(end[0]-start[0]);
        double py=start[1]+tp*(end[1]-start[1]);
        return dis_points(xp,yp,px,py);
    }
}

double* actin::get_intpoint(double xp, double yp)
{
    double* coordinates = new double[2];
    double l2 = pow(dis_points(start[0], start[1], end[0], end[1]) , 2);
    if (l2==0) {
        coordinates[0]=start[0];
        coordinates[1]=start[1];
    }
    double tp=dot(xp-start[0],yp-start[1],end[0]-start[0],end[1]-start[1])/l2;
    if (tp<0) {
        coordinates[0]=start[0];
        coordinates[1]=start[1];
    }
    else if(tp>1.0){
        coordinates[0]=end[0];
        coordinates[1]=end[1];
    }
    else{
        coordinates[0]=start[0]+tp*(end[0]-start[0]);
        coordinates[1]=start[1]+tp*(end[1]-start[1]);
    }
    return coordinates;
}

double actin::get_int_angle(double xp, double yp)
{
    double angle;
    double xcor,ycor;
    double slope=(end[1]-start[1])/(end[0]-start[0]);
    double yintercept=y-slope*x;
    xcor=(slope*yp + xp - slope*yintercept)/(slope*slope + 1);
    ycor=(slope*slope*yp + slope*xp + yintercept)/(1 + slope*slope);
    angle=atan2((ycor-yp),(xcor-xp));
    return angle;
}

double*  actin::get_direction()
{
    return e;
}

double actin::get_length()
{
    return ld;
}

double* actin::get_forces()
{
    return forces;
}

void actin::update_force(double f1, double f2, double f3)
{
    forces[0]+=f1;
    forces[1]+=f2;
    forces[2]+=f3;
}

double* actin::get_friction()
{
    double* fric = new double[3];
    fric[0]=2*pi*a_vis*ld/log(ld/diameter);
    fric[1]=2*fric[0];
    fric[2]=fric[0]*pow(ld,2)/6;
    return fric;
}

double actin::get_xcm()
{
    return x;
}

double actin::get_ycm()
{
    return y;
}

double actin::get_angle()
{
    return phi;
}

double * actin::get_start(){
    return start;
}

double * actin::get_end(){
    return end;
}

void actin::set_xcm(double xcm)
{
    x = xcm;
}

void actin::set_ycm(double ycm)
{
    y = ycm;
}

void actin::set_phi(double theta)
{
    phi = theta;
}

std::vector<std::vector<int> > actin::get_quadrants()
{ 
    return quad; 
}
