/*
 * motor.cpp
 *  
 *
 *  Created by Shiladitya Banerjee on 9/3/13.
 *  Copyright 2013 University of Chicago. All rights reserved.
 *
 */

#include "globals.h"
#include "motor.h"
#include "actin_ensemble.h"
//#include "actin.h"

//motor class
motor::motor(double mx, double my, double mang, double mlen, actin_ensemble* network, int state0, int state1, int
        aindex0, int aindex1, double fovx, double fovy, double delta_t, double temp,
        double v0, double stiffness, double ron, double roff, double
        rend, double actin_len, double vis, std::string col) {
    vs=v0;//rng_n(v0,0.4);//rng(v0-0.3,v0+0.3);
    dm=0.25;//actin_len/10;
    mk=stiffness;//rng(10,100);
    fmax=mk*dm*2;//rng(1,20);
    mld=mlen;
    kon=ron;
    koff=roff;
    kend=rend;
    mphi=mang;
    dt = delta_t;
    temperature = temp;
    hx[0]=mx-0.5*mld*cos(mphi);
    hy[0]=my-0.5*mld*sin(mphi);
    hx[1]=hx[0]+mld*cos(mphi);
    hy[1]=hy[0]+mld*sin(mphi);
    mobility=log(10)/(4*pi*vis*mld);
    state[0]=state0;
    state[1]=state1;
    aindex[0]=aindex0;//   actin index for head in state[0]
    aindex[1]=aindex1;// actin index for head in state[1]
    actin_network=network;
    //		pos_actin[0]=0; // I don't think this variable get's ACCESSED anywhere 
    pos_a_end[0]=0; //distance from pointy end -- by default 0
    pos_a_end[1]=0;

    if (state0){
        pos_a_end[0] = fmin(dis_points(hx[0],hy[0],actin_network->get_start(aindex[0])[0],actin_network->get_start(aindex[0])[1]),
                dis_points(hx[0],hy[0],actin_network->get_end(aindex[0])[0],actin_network->get_end(aindex[0])[1]));
    }
    if (state1){
        pos_a_end[1] = fmin(dis_points(hx[1],hy[1],actin_network->get_start(aindex[1])[0],actin_network->get_start(aindex[1])[1]),
                dis_points(hx[1],hy[1],actin_network->get_end(aindex[1])[0],actin_network->get_end(aindex[1])[1]));
    }
    
    fov[0]=fovx;
    fov[1]=fovy;
    color = col; 

}

motor::~motor(){ 
};

//return motor state with a given head number
int* motor::get_states() 
{
    return state;
}

double* motor::get_hx()
{
    return hx;
}

double* motor::get_hy()
{
    return hy;
}

std::string motor::get_color()
{
    return color;
}

double motor::tension()
{
    double lf=dis_points(hx[0],hy[0],hx[1],hy[1]);
    return mk*(lf-mld)/fmax;
}

//check for attachment of unbound heads given head index (0 for head 1, and 1 for head 2)
void motor::attach(int hd)
{
    dist.clear();
    dist=actin_network->get_dist(hx[hd],hy[hd]);
    if(!dist.empty()){
        for (std::map<int,double>::iterator it=dist.begin(); it!=dist.end(); ++it)
        { 
            if (it->second <= dm && aindex[pr(hd)]!=it->first) {
                onrate=kon*exp(-((it->second)*(it->second))/(dm*dm));
                if (event(onrate,dt)==1) {
                    //update state
                    state[hd]=1;
                    aindex[hd]=it->first;
                    //update head position
                    
                    double * intpoint = actin_network->get_intpoints(it->first,hx[hd],hy[hd]);
                    hx[hd] = intpoint[0];
                    hy[hd] = intpoint[1];
                    delete intpoint;

                    if (state[pr(hd)]==0) {
                        hx[pr(hd)]=hx[hd]+pow(-1,hd)*mld*cos(mphi);
                        hy[pr(hd)]=hy[hd]+pow(-1,hd)*mld*sin(mphi);
                    }
                    else {
                        mphi=atan2((hy[1]-hy[0]),(hx[1]-hx[0]));
                    }

                    //                        pos_actin[hd]=dis_points(hx[hd],hy[hd],actin_network->get_xcm(aindex[hd]),actin_network->get_ycm(aindex[hd]));
                    pos_a_end[hd]=dis_points(hx[hd],hy[hd],actin_network->get_end(aindex[hd])[0],actin_network->get_end(aindex[hd])[1]);
                    break;
                }
            }
        }
    }	
} 

//perform brownian motion and shear if head unattached
void motor::brownian(double t, double gamma)
{
    if (state[0]==0 && state[1]==0) {

        xm[0]=hx[0]+sqrt(dt*mobility*2*temperature)*rng_n(0,1) - mk*dt*(hx[0]-hx[1]+mld*cos(mphi)) + gamma*dt*hy[0];
        xm[1]=hx[1]+sqrt(dt*mobility*2*temperature)*rng_n(0,1) + mk*dt*(hx[0]-hx[1]+mld*cos(mphi)) + gamma*dt*hy[1];
        ym[0]=hy[0]+sqrt(dt*mobility*2*temperature)*rng_n(0,1) - mk*dt*(hy[0]-hy[1]+mld*sin(mphi));
        ym[1]=hy[1]+sqrt(dt*mobility*2*temperature)*rng_n(0,1) + mk*dt*(hy[0]-hy[1]+mld*sin(mphi));
        reflect(t, gamma, xm[0],xm[1],ym[0],ym[1]);
        mphi=atan2((hy[1]-hy[0]),(hx[1]-hx[0]));
    }
    else if (state[0]==0 || state[1]==0) {
        int hd=state[0];//magically equivalent to 
                        //int hd = (hd for which state[hd] = 0)

        xm[hd]=hx[hd]+sqrt(dt*mobility*2*temperature)*rng_n(0,1) - mk*dt*(hx[hd]-hx[pr(hd)]+pow(-1,hd)*mld*cos(mphi)) + gamma*dt*hy[hd];
        ym[hd]=hy[hd]+sqrt(dt*mobility*2*temperature)*rng_n(0,1) - mk*dt*(hy[hd]-hy[pr(hd)]+pow(-1,hd)*mld*sin(mphi));
        xm[pr(hd)]=hx[pr(hd)];
        ym[pr(hd)]=hy[pr(hd)];
        reflect(t, gamma, xm[0],xm[1],ym[0],ym[1]);
        mphi=atan2((hy[1]-hy[0]),(hx[1]-hx[0]));
    }
    else {
        return;
    }

}

//stepping and detachment kinetics of a single bound head 
void motor::step_onehead(int hd)
{

    if (event(koff,dt)==1) {
        state[hd]=0;
        aindex[hd]=-1;
        //			pos_actin[hd]=0;
        pos_a_end[hd]=0;
        hx[hd]=hx[pr(hd)]-pow(-1,hd)*mld*cos(mphi);
        hy[hd]=hy[pr(hd)]-pow(-1,hd)*mld*sin(mphi);

    }
    else
    {
        pos_temp=pos_a_end[hd]+dt*vs;
        move_end_detach(hd,pos_temp);


    }
}

//stepping and detachment kinetics of a motor with both heads attached
void motor::step_twoheads()
{
    stretch=dis_points(hx[0],hy[0],hx[1],hy[1])-mld;
    //fm = vec(Fm).(-vec(u)) 
    fm[0]=mk*((hx[0]-hx[1]+mld*cos(mphi))*actin_network->get_direction(aindex[0])[0] +
            (hy[0]-hy[1]+mld*sin(mphi))*actin_network->get_direction(aindex[0])[1]); 
    vm[0]=velocity(vs,fm[0],fmax);
    fm[1]=mk*(-(hx[0]-hx[1]+mld*cos(mphi))*actin_network->get_direction(aindex[1])[0] -
            (hy[0]-hy[1]+mld*sin(mphi))*actin_network->get_direction(aindex[1])[1]);

    vm[1]=velocity(vs,fm[1],fmax);
    /*impose force-dependent bell's law on detachment rates*/
    offrate[0]=koff*exp(fabs(fm[0])/fmax);
    offrate[1]=koff*exp(fabs(fm[1])/fmax);

    if (event(offrate[0],dt)==1) {
        state[0]=0;
        aindex[0]=-1;
        hx[0]=hx[1]-mld*cos(mphi);
        hy[0]=hy[1]-mld*sin(mphi);
        //			pos_actin[0]=0;
        pos_a_end[0]=0;
        pos_temp=pos_a_end[1]+dt*vs;
        move_end_detach(1,pos_temp);
    }
    else if (event(offrate[1],dt)==1) {
        state[1]=0;
        aindex[1]=-1;
        hx[1]=hx[0]+mld*cos(mphi);
        hy[1]=hy[0]+mld*sin(mphi);
        //            pos_actin[1]=0;
        pos_a_end[1]=0;
        pos_temp=pos_a_end[0]+dt*vs;
        move_end_detach(0,pos_temp);
    }

    else {

        pos_temp=pos_a_end[0]+dt*vm[0];
        move_end_detach(0,pos_temp);
        pos_temp=pos_a_end[1]+dt*vm[1];
        move_end_detach(1,pos_temp); 
    }
}


void motor::actin_update()
{
    if (state[0]==1 && state[1]==1) {
        stretch=dis_points(hx[0],hy[0],hx[1],hy[1])-mld;
//        std::cout<<"DEBUG: actin_update: color = "<< color<< "\tstretch = "<<stretch<<"\n";
        forcex[0]=-mk*(hx[0]-hx[1]+mld*cos(mphi));
        forcex[1]=-forcex[0];
        forcey[0]=-mk*(hy[0]-hy[1]+mld*sin(mphi));
        forcey[1]=-forcey[0];
        force_par[0]=forcex[0]*actin_network->get_direction(aindex[0])[0] + forcey[0]*actin_network->get_direction(aindex[0])[1];
        force_perp[0]=-forcex[0]*actin_network->get_direction(aindex[0])[1] + forcey[0]*actin_network->get_direction(aindex[0])[0];
        force_par[1]=forcex[1]*actin_network->get_direction(aindex[1])[0] + forcey[1]*actin_network->get_direction(aindex[1])[1];
        force_perp[1]=-forcex[1]*actin_network->get_direction(aindex[1])[1] + forcey[1]*actin_network->get_direction(aindex[1])[0];

        torque[0]=cross(hx[0]-actin_network->get_xcm(aindex[0]),hy[0]-actin_network->get_ycm(aindex[0]),forcex[0],forcey[0]);
        torque[1]=cross(hx[1]-actin_network->get_xcm(aindex[1]),hy[1]-actin_network->get_ycm(aindex[1]),forcex[1],forcey[1]);
        actin_network->update_forces(aindex[0],force_par[0],force_perp[0],torque[0]);
        actin_network->update_forces(aindex[1],force_par[1],force_perp[1],torque[1]);
    }
    else
        return;

}

void motor::update_shape()
{
    if (state[0]==1 && state[1]==1) {
        hx[0]=actin_network->get_end(aindex[0])[0]-pos_a_end[0]*actin_network->get_direction(aindex[0])[0];
        hy[0]=actin_network->get_end(aindex[0])[1]-pos_a_end[0]*actin_network->get_direction(aindex[0])[1];
        //            pos_actin[0]=dis_points(hx[0],hy[0],actin_network->get_xcm(aindex[0]),actin_network->get_ycm(aindex[0]));
        hx[1]=actin_network->get_end(aindex[1])[0]-pos_a_end[1]*actin_network->get_direction(aindex[1])[0];
        hy[1]=actin_network->get_end(aindex[1])[1]-pos_a_end[1]*actin_network->get_direction(aindex[1])[1];
        //            pos_actin[1]=dis_points(hx[1],hy[1],actin_network->get_xcm(aindex[1]),actin_network->get_ycm(aindex[1]));
        mphi=atan2((hy[1]-hy[0]),(hx[1]-hx[0]));
    }
    else if(state[0]==1 && state[1]==0)
    {
        hx[0]=actin_network->get_end(aindex[0])[0]-pos_a_end[0]*actin_network->get_direction(aindex[0])[0];
        hy[0]=actin_network->get_end(aindex[0])[1]-pos_a_end[0]*actin_network->get_direction(aindex[0])[1];
        //            pos_actin[0]=dis_points(hx[0],hy[0],actin_network->get_xcm(aindex[0]),actin_network->get_ycm(aindex[0]));
        mphi=atan2((hy[1]-hy[0]),(hx[1]-hx[0]));
    }
    else if(state[0]==0 && state[1]==1)
    {
        hx[1]=actin_network->get_end(aindex[1])[0]-pos_a_end[1]*actin_network->get_direction(aindex[1])[0];
        hy[1]=actin_network->get_end(aindex[1])[1]-pos_a_end[1]*actin_network->get_direction(aindex[1])[1];
        //            pos_actin[1]=dis_points(hx[1],hy[1],actin_network->get_xcm(aindex[1]),actin_network->get_ycm(aindex[1]));
        mphi=atan2((hy[1]-hy[0]),(hx[1]-hx[0]));
    }
    else
    {
        return;
    }


}

inline void motor::move_end_detach(int hd, double pos)
{
    double rod_length = actin_network->get_alength(aindex[hd]);
    if (pos >= rod_length) { 
        
        if (actin_network->is_polymer_start(aindex[hd])){
            if (event(kend,dt)==1) {
                state[hd]=0;
                aindex[hd]=-1;

                pos_a_end[hd]=0;
                hx[hd]=hx[pr(hd)]-pow(-1,hd)*mld*cos(mphi);
                hy[hd]=hy[pr(hd)]-pow(-1,hd)*mld*sin(mphi);
            }
            else {
                hx[hd]=actin_network->get_end(aindex[hd])[0]-pos_a_end[hd]*actin_network->get_direction(aindex[hd])[0];
                hy[hd]=actin_network->get_end(aindex[hd])[1]-pos_a_end[hd]*actin_network->get_direction(aindex[hd])[1];
                mphi=atan2((hy[1]-hy[0]),(hx[1]-hx[0]));
            }
        }
        else{ 
            /*Move the motor to the next rod on the filament
             *At the projected new position along that filament
             */
            
            double rod_angle = actin_network->get_angle(aindex[hd]);

            aindex[hd] = aindex[hd] - 1;
            double new_rod_angle = actin_network->get_angle(aindex[hd]);
            
            pos_a_end[hd] = (pos - rod_length)*cos(new_rod_angle - rod_angle);
            hx[hd]=actin_network->get_end(aindex[hd])[0]-pos_a_end[hd]*actin_network->get_direction(aindex[hd])[0];
            hy[hd]=actin_network->get_end(aindex[hd])[1]-pos_a_end[hd]*actin_network->get_direction(aindex[hd])[1];
            mphi=atan2((hy[1]-hy[0]),(hx[1]-hx[0]));
    
        }
    }
    else {
        pos_a_end[hd]=pos;
        hx[hd]=actin_network->get_end(aindex[hd])[0]-pos_a_end[hd]*actin_network->get_direction(aindex[hd])[0];
        hy[hd]=actin_network->get_end(aindex[hd])[1]-pos_a_end[hd]*actin_network->get_direction(aindex[hd])[1];
        mphi=atan2((hy[1]-hy[0]),(hx[1]-hx[0]));
    }
    
    stretch=dis_points(hx[0],hy[0],hx[1],hy[1])-mld;

}

inline void motor::reflect(double t, double gamma, double x1, double x2, double y1, double y2)
{
    //Calculate the sheared simulation bounds (at this height)
    double xleft, xright;
    xleft  = std::max(-fov[0] * 0.5 + gamma * y1 * t, -fov[0] * 0.5 + gamma * y2 * t);
    xright = std::min( fov[0] * 0.5 + gamma * y1 * t,  fov[0] * 0.5 + gamma * y2 * t);
    if (xleft < x1 && x1 < xright
            && xleft < x2 && x2 < xright
            && -fov[1]*0.5 < y1 && y1 < fov[1]*0.5 
            && -fov[1]*0.5 < y2 && y2 < fov[1]*0.5) {
        hx[0]=x1;
        hx[1]=x2;
        hy[0]=y1;
        hy[1]=y2;
    }
    else if (x1>=xright || x1<=xleft)
    {
        hx[1]=x2;
        hy[0]=y1;
        hy[1]=y2;   
    }
    else if (x2>=xright || x2<=xleft)
    {
        hx[0]=x1;
        hy[0]=y1;
        hy[1]=y2;
    }
    else if(y1>=fov[1]*0.5 || y1<=fov[1]*0.5)
    {
        hx[0]=x1;
        hx[1]=x2;
        hy[1]=y2;
    }
    else{
        hx[0]=x1;
        hx[1]=x2;
        hy[0]=y1;
    }
}
