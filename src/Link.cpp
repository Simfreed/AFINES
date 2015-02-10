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

Link::Link(double len, double stretching_stiffness, double bending_stiffness, 
        filament* f, int aindex0, int aindex1)
{
    kl              =   stretching_stiffness;
    kb              =   bending_stiffness;
    ld              =   len;
    aindex[0]       =   aindex0;
    aindex[1]       =   aindex1;
    fil             =   f;

    // initialize the coordinates of the heads:
    hx[0] = 0;
    hx[1] = 0;
    hy[0] = 0;
    hy[1] = 0;

    this->step();
}
Link::~Link(){ 
    //std::cout<<"DELETING LINK\n";
};

double* Link::get_hx(){
    return hx;
}

double* Link::get_hy(){
    return hy;
}

void Link::set_aindex1(int i){
    aindex[1] = i;
    this->step();
}

// stepping kinetics
void Link::step()
{
    //CONVENTION: head 0 will be connected to the POINTY end of a filament
    //            head 1 will be connected to the BARBED end of a filament
    if (aindex[0]==-1){ //leftmost end of the polymer
        double * start1 = fil->get_rod(aindex[1])->get_start();
        hx[1] = start1[0];
        hy[1] = start1[1];
        hx[0] = hx[1] - ld*cos( fil->get_rod(aindex[1])->get_angle() );
        hy[0] = hy[1] - ld*sin( fil->get_rod(aindex[1])->get_angle() );
    }else if(aindex[1] == -1){ //rightmost end of the polymer
        double * end0 = fil->get_rod(aindex[0])->get_end();
        
        hx[0] = end0[0];
        hy[0] = end0[1];
        hx[1] = hx[0] + ld*cos( fil->get_rod(aindex[0])->get_angle() );
        hy[1] = hy[0] + ld*sin( fil->get_rod(aindex[0])->get_angle() );
    }else{
        double * end0 = fil->get_rod(aindex[0])->get_end();
        double * start1 = fil->get_rod(aindex[1])->get_start();
         
        hx[0] = end0[0];
        hy[0] = end0[1];
        hx[1] = start1[0]; 
        hy[1] = start1[1];
    }

    xcm = (hx[0]+hx[1])/2.0;
    ycm = (hy[0]+hy[1])/2.0;
    phi=atan2(hy[1]-hy[0],hx[1]-hx[0]);

}

double Link::get_stretch_force(){
    //cout<<"\nDEBUG:Stretch force = "<<kl*(dis_points(hx[0],hy[0],hx[1],hy[1])-ld);
    return kl * (dis_points(hx[0],hy[0],hx[1],hy[1])-ld);
}

void Link::filament_update()
{

    double force_stretch = this->get_stretch_force();
    double * e0, * e1;
    
    //cout<<"\nDEBUG:stretch = "<<force_stretch;
    if (aindex[0] != -1){
        e0 = fil->get_rod(aindex[0])->get_direction();
        forcex[0]       =  -force_stretch * cos(phi); 
        forcey[0]       =  -force_stretch * sin(phi); 
        force_par[0]    =   forcex[0]*e0[0] + forcey[0]*e0[1];
        force_perp[0]   =  -forcex[0]*e0[1] + forcey[0]*e0[0];
        torque[0]       =   cross(hx[0]-fil->get_rod(aindex[0])->get_xcm(),hy[0]-fil->get_rod(aindex[0])->get_ycm(),forcex[0],forcey[0]);
        //cout<<"\nDEBUG: rod "<<aindex[0]<<" : ("<<force_par[0]<<","<<force_perp[0]<<","<<torque[0]<<")";
        fil->update_forces(aindex[0],force_par[0],force_perp[0],torque[0]);
    }

    if (aindex[1] != -1){
        e1 = fil->get_rod(aindex[1])->get_direction();
        forcex[1]       =    force_stretch * cos(phi);
        forcey[1]       =    force_stretch * sin(phi);
        force_par[1]    =    forcex[1]*e1[0] + forcey[1]*e1[1];
        force_perp[1]   =   -forcex[1]*e1[1] + forcey[1]*e1[0];
        torque[1]       =   cross(hx[1]-fil->get_rod(aindex[1])->get_xcm(),hy[1]-fil->get_rod(aindex[1])->get_ycm(),forcex[1],forcey[1]);
        //cout<<"\nDEBUG: rod "<<aindex[1]<<" : ("<<force_par[1]<<","<<force_perp[1]<<","<<torque[1]<<")";
        fil->update_forces(aindex[1],force_par[1],force_perp[1],torque[1]);
    }

}

/*********************************************************
 * Calculates the x and y forces on the end of rod due   *
 *      to the link that connects them                   *
 * For Link j in the filament, this function             *
 *      returns an array of doubles callled forces       *
 * where:                                                *
 *      forces[0] is the x force on the end of rod j-1   *
 *      forces[1] is the y force on the end of rod j-1   *
 *      forces[2] is the x force at the start of rod j   *
 *      forces[3] is the y force at the start of rod j   *
 *********************************************************/

double* Link::get_forces()
{

    double* forces = new double[4]; 
    
    if (aindex[0] == -1 || aindex[1] == -1){
        forces[0] = 0;
        forces[1] = 0;
        forces[2] = 0;
        forces[3] = 0;
    }
    else
    {
        double force_stretch = this->get_stretch_force();
    //    cout<<"\nDEBUG: force_stretch = "<<force_stretch; 
        forces[0] = force_stretch * cos(phi); 
        forces[1] = force_stretch * sin(phi); 
        forces[2] = -force_stretch * cos(phi);
        forces[3] = -force_stretch * sin(phi);
    /*    cout<<"\nDEBUG: forces[0] = "<<forces[0]; 
        cout<<"\nDEBUG: forces[1] = "<<forces[1]; 
        cout<<"\nDEBUG: forces[2] = "<<forces[2]; 
        cout<<"\nDEBUG: forces[3] = "<<forces[3]; 
    */
    }
    
    return forces;
}

double Link::get_kb(){
    return kb;
}

double Link::get_kl(){
    return kl;
}

double Link::get_length(){
    return ld;
}
double Link::get_xcm(){
    return xcm;
}

double Link::get_ycm(){
    return ycm;
}

std::string Link::write(){
    return std::to_string(hx[0]) + "\t" + std::to_string(hy[0]) + "\t" + std::to_string(hx[1]-hx[0]) + "\t" 
        + std::to_string(hy[1]-hy[0]) + "\n";
}

std::string Link::to_string(){
    
    char buffer [100];
    sprintf(buffer, "aindex[0] = %d;\t aindex[1] = %d;\t kl = %f;\t kb = %f;\t ld = %f\nfilament : \n",
                        aindex[0], aindex[1], kl, kb, ld);

    return buffer + fil->to_string();

}

bool Link::operator==(const Link& that) 
{
    /*Note: you can't compare the filament objects because that will lead to infinite recursion;
     * this function requires the filament poiner to be identical to evaluate to true*/
    return (this->aindex[0] == that.aindex[0] && this->aindex[1] == that.aindex[1] &&
            this->kl == that.kl && this->kb == that.kb &&
            this->ld == that.ld && this->fil == that.fil);
}

bool Link::is_similar(const Link& that) 
{
    
    /* Same as ==; but doesn't compare the filament pointer*/

    return (this->aindex[0] == that.aindex[0] && this->aindex[1] == that.aindex[1] &&
            this->kl == that.kl && this->kb == that.kb &&
            this->ld == that.ld);
}


MidLink::MidLink(double len, double stretching_stiffness, double bending_stiffness, 
        filament* f, int aindex0, int aindex1)
    : Link(len, stretching_stiffness, bending_stiffness, f, aindex0, aindex1)
{
    // initialize the coordinates of the heads:
    is_linear=false;
    this->step();
}
void MidLink::step()
{

    hx[0]=fil->get_rod(aindex[0])->get_xcm();
    hy[0]=fil->get_rod(aindex[0])->get_ycm();
    hx[1]=fil->get_rod(aindex[1])->get_xcm();
    hy[1]=fil->get_rod(aindex[1])->get_ycm();

    phi=atan2(hy[1]-hy[0],hx[1]-hx[0]);
    
    xcm = (hx[0]+hx[1])/2.0;
    ycm = (hy[0]+hy[1])/2.0;

}

double MidLink::get_stretch_force(){
    double stretch = dis_points(hx[0],hy[0],hx[1],hy[1])-ld; 
    if (is_linear){
        return kb * stretch;
    }
    else{
        if (stretch > 0)
            return kb;
        else if(stretch < 0)
            return -1*kb;
        else 
            return 0;
    }
}

void MidLink::set_linear(bool linear){
    is_linear = linear;
}
