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
actin::actin(){}

actin::actin(double xcm, double ycm, double angle, double len, double fovx, double fovy, int nx, int ny, double vis)
{
    x=xcm;
    y=ycm;
    phi=angle;
    ld=len;
    diameter= 0.006; //source: Google //ld/40;
    a_vis=vis;
    
    fov[0] = fovx;
    fov[1] = fovy;
    nq[0]  = nx;
    nq[1]  = ny;
    
    this->update();
}

actin::actin(const actin& other){
    
    x = other.x;
    y = other.y;
    phi = other.phi;
    ld = other.ld;
    diameter = other.diameter;
    a_vis = other.a_vis;
    fov[0] = other.fov[0];
    fov[1] = other.fov[1];
    nq[0] = other.nq[0];
    nq[1] = other.nq[1];

    this->update();
}

actin::~actin(){ 
    //std::cout<<"DELETING ACTIN\n";
};

// Updates all derived quantities of a monomer
void actin::update(){
    
    start[0]=x-ld*0.5*cos(phi);
    start[1]=y-ld*0.5*sin(phi);
    end[0]=x+ld*0.5*cos(phi);
    end[1]=y+ld*0.5*sin(phi);
    
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
        lower_limit = int(floor(start[0]/fov[0]*nq[0]));
        if(lower_limit > 0){lower_limit--;};
        upper_limit = int(ceil(end[0]/fov[0]*nq[0]));
        if(upper_limit < nq[0]-1){upper_limit++;};

        for(index = lower_limit; index < upper_limit;index++){tmp.push_back(index);};
    }
    else
    {
        lower_limit = int(floor(end[0]/fov[0]*nq[0]));
        if(lower_limit > 0){lower_limit--;};
        upper_limit = int(ceil(start[0]/fov[0]*nq[0]));
        if(upper_limit < nq[0]-1){upper_limit++;};
        for(index = lower_limit; index < upper_limit;index++){tmp.push_back(index);};
    };
    quad.push_back(tmp);

    //quadrant numbers crossed by the actin in y-direction
    tmp.clear();
    if(start[1] <= end[1])
    {
        lower_limit = int(floor(start[1]/fov[1]*nq[1]));
        if(lower_limit > 0){lower_limit--;};
        upper_limit = int(ceil(end[1]/fov[1]*nq[1]));
        if(upper_limit < nq[1]-1){upper_limit++;};
        for(index = lower_limit; index < upper_limit;index++){tmp.push_back(index);}
    }
    else
    {
        lower_limit = int(floor(end[1]/fov[1]*nq[1]));
        if(lower_limit > 0){lower_limit--;};
        upper_limit = int(ceil(start[1]/fov[1]*nq[1]));
        if(upper_limit < nq[1]-1){upper_limit++;};
        for(index = lower_limit; index < upper_limit;index++){tmp.push_back(index);}
    };
    quad.push_back(tmp);

}
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
    fric[2]=fric[0]*pow(ld,2)/4;
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

bool actin::operator==(const actin& that) 
{
    double err = eps; 
    return (close(x , that.x , err) && close(y , that.y , err) &&
            close(phi , that.phi , err) && close(ld , that.ld , err) &&
            close(a_vis , that.a_vis , err) && close(forces[0] , that.forces[0] , err) &&
            close(forces[1] , that.forces[1] , err) && close(forces[2] , that.forces[2], err) 
           );
}

std::string actin::write()
{
    return std::to_string(start[0]) + "\t" + std::to_string(start[1]) + "\t" + 
           std::to_string(end[0]-start[0]) + "\t" + std::to_string(end[1]-start[1]) + "\n";
 
}

std::string actin::to_string()
{
    return "x : " + std::to_string(x) + "\ty : " + std::to_string(y) + "\tphi : " + 
           std::to_string(phi) + "\tld : " + std::to_string(ld) + "\ta_vis : " +
           std::to_string(a_vis) + "\tforces[0] : " + std::to_string(forces[0]) + "\tforces[1] : " +
           std::to_string(forces[1]) + "\tforces[2] : " + std::to_string(forces[2]) + "\n";
 
}

// EXCLUDED VOLUME STUFF
//

void actin::set_gay_berne(double s0, double e0, double epsS, double epsE, double m, double n)
{
    sigma0 = s0;
    eps0 = e0;
    mu = m;
    nu = n;
    
    chi      = (ld * ld - diameter * diameter)  / (ld * ld + diameter * diameter);
    chiPrime = (pow(epsS, 1.0/mu)-pow(epsE, 1.0/mu))/(pow(epsS, 1.0/mu) + pow(epsE, 1.0/mu));

}

double actin::get_sigma0()
{
    return sigma0;
}

double actin::get_eps0()
{
    return eps0;
}

double actin::get_chi()
{
    return chi;
}

double actin::get_chiPrime()
{
    return chiPrime;
}

//Calculates the Gay Berne forces and torques between two actin rods

double * actin::calc_gay_berne(actin * a)
{
    double u1u2, ru1, ru2, rx, ry, R, sigma, eps, epsPrime, epsUnPrime;
    double * u1, * u2;
    double * forces = new double[3]; 
    u1 = e;
    u2 = a->get_direction();
    rx = x - a->get_xcm();
    ry = y - a->get_ycm();
    R = rx*rx + ry*ry;
    
    u1u2 = dot(u1[0], u1[1], u2[0], u2[1]);
    ru1 = dot(rx, ry, u1[0], u1[1]);
    ru2 = dot(rx, ry, u1[0], u1[1]);
    
    epsUnPrime = eps0 * pow((1.0-pow(chi,2)*u1u2*u1u2),-0.5);
    epsPrime   = 1 - chiPrime / 2.0 * ( pow(ru1 + ru2, 2)/(1+chiPrime*u1u2) + pow(ru1 - ru2, 2)/(1-chiPrime*u1u2) );

    return forces;

}

double * actin::get_fov(){
    return fov;
}

double * actin::get_nq(){
    return nq;
}

double actin::get_viscosity(){
    return a_vis;
}

double actin::get_diameter(){
    return diameter;
}
