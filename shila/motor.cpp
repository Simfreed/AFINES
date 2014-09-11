/*
 * motor.cpp
 *  
 *
 *  Created by Shiladitya Banerjee on 9/3/13.
 *  Copyright 2013 University of Chicago. All rights reserved.
 *
 */

#include "generic_functions.cpp"
#include "motor.h"
#include "actin_ensemble.h"
#include "actin.h"

//motor class
motor::motor(double mx, double my, double mang, double mlen, actin_ensemble* network, int state0, int state1, int
        aindex0, int aindex1, double fovx, double fovy, double v0, double stiffness, double ron, double roff, double
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
        pos_a_end[0] = fmin(dis_points(hx[0],hy[0],actin_network->get_ends(aindex[0])[0],actin_network->get_ends(aindex[0])[1]),
                dis_points(hx[0],hy[0],actin_network->get_ends(aindex[0])[2],actin_network->get_ends(aindex[0])[3]));
        std::cout<<"DEBUG: distance of head 0 from pointy end of actin monomer "<<aindex[0]<<" : "<<pos_a_end[0]<<"\n";
    }
    if (state1){
        pos_a_end[1] = fmin(dis_points(hx[1],hy[1],actin_network->get_ends(aindex[1])[0],actin_network->get_ends(aindex[1])[1]),
                dis_points(hx[1],hy[1],actin_network->get_ends(aindex[1])[2],actin_network->get_ends(aindex[1])[3]));
        std::cout<<"DEBUG: distance of head 1 from pointy end of actin monomer "<<aindex[1]<<" : "<<pos_a_end[0]<<"\n";
    }
    //        pos_actin[1]=0;
    fov[0]=fovx;
    fov[1]=fovy;
    color = col; 

}

//return motor state with a given head number
int* motor::get_states() 
{
    int* sptr;
    sptr=state;
    return sptr;
}


double* motor::get_heads()
{
    double h[4];
    double *gh;
    h[0]=hx[0];
    h[1]=hy[0];
    h[2]=hx[1];
    h[3]=hy[1];
    gh=h;
    return gh;
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
                    hx[hd]=actin_network->get_intpoints(it->first,hx[hd],hy[hd])[0];
                    hy[hd]=actin_network->get_intpoints(it->first,hx[hd],hy[hd])[1];
                    if (state[pr(hd)]==0) {
                        hx[pr(hd)]=hx[hd]+pow(-1,hd)*mld*cos(mphi);
                        hy[pr(hd)]=hy[hd]+pow(-1,hd)*mld*sin(mphi);
                    }
                    else {
                        mphi=atan2((hy[1]-hy[0]),(hx[1]-hx[0]));
                    }

                    //                        pos_actin[hd]=dis_points(hx[hd],hy[hd],actin_network->get_position(aindex[hd])[0],actin_network->get_position(aindex[hd])[1]);
                    pos_a_end[hd]=dis_points(hx[hd],hy[hd],actin_network->get_ends(aindex[hd])[2],actin_network->get_ends(aindex[hd])[3]);
                    break;
                }
            }
        }
    }	
} 

//perform brownian motion if head unattached
void motor::brownian()
{
    if (state[0]==0 && state[1]==0) {

        xm[0]=hx[0]+sqrt(dt*mobility*6*temperature)*rng_n(0,1) - mk*dt*(hx[0]-hx[1]+mld*cos(mphi));
        xm[1]=hx[1]+sqrt(dt*mobility*6*temperature)*rng_n(0,1) + mk*dt*(hx[0]-hx[1]+mld*cos(mphi));
        ym[0]=hy[0]+sqrt(dt*mobility*6*temperature)*rng_n(0,1) - mk*dt*(hy[0]-hy[1]+mld*sin(mphi));
        ym[1]=hy[1]+sqrt(dt*mobility*6*temperature)*rng_n(0,1) + mk*dt*(hy[0]-hy[1]+mld*sin(mphi));
        reflect(xm[0],xm[1],ym[0],ym[1]);
        mphi=atan2((hy[1]-hy[0]),(hx[1]-hx[0]));
    }
    else if (state[0]==0 || state[1]==0) {
        int hd=state[0];

        xm[hd]=hx[hd]+sqrt(dt*mobility*6*temperature)*rng_n(0,1) - mk*dt*(hx[hd]-hx[pr(hd)]+pow(-1,hd)*mld*cos(mphi));
        ym[hd]=hy[hd]+sqrt(dt*mobility*6*temperature)*rng_n(0,1) - mk*dt*(hy[hd]-hy[pr(hd)]+pow(-1,hd)*mld*sin(mphi));
        xm[pr(hd)]=hx[pr(hd)];
        ym[pr(hd)]=hy[pr(hd)];
        reflect(xm[0],xm[1],ym[0],ym[1]);
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
        move_end_detach(hd,vs,pos_temp);


    }
}

//stepping and detachment kinetics of a motor with both heads attached
void motor::step_twoheads()
{
    stretch=dis_points(hx[0],hy[0],hx[1],hy[1])-mld;
    std::cout<<"DEBUG: step_twoheads: color = "<< color<< "\tstretch = "<<stretch<<"\n";
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
        move_end_detach(1,vs,pos_temp);
    }
    else if (event(offrate[1],dt)==1) {
        state[1]=0;
        aindex[1]=-1;
        hx[1]=hx[0]+mld*cos(mphi);
        hy[1]=hy[0]+mld*sin(mphi);
        //            pos_actin[1]=0;
        pos_a_end[1]=0;
        pos_temp=pos_a_end[0]+dt*vs;
        move_end_detach(0,vs,pos_temp);
    }

    else {

        pos_temp=pos_a_end[0]+dt*vm[0];
        move_end_detach(0,vm[0],pos_temp);
        pos_temp=pos_a_end[1]+dt*vm[1];
        move_end_detach(1,vm[1],pos_temp); 
    }
}


void motor::actin_update()
{
    if (state[0]==1 && state[1]==1) {
        stretch=dis_points(hx[0],hy[0],hx[1],hy[1])-mld;
        std::cout<<"DEBUG: actin_update: color = "<< color<< "\tstretch = "<<stretch<<"\n";
        forcex[0]=-mk*(hx[0]-hx[1]+mld*cos(mphi));
        forcex[1]=-forcex[0];
        forcey[0]=-mk*(hy[0]-hy[1]+mld*sin(mphi));
        forcey[1]=-forcey[0];
        force_par[0]=forcex[0]*actin_network->get_direction(aindex[0])[0] + forcey[0]*actin_network->get_direction(aindex[0])[1];
        force_perp[0]=-forcex[0]*actin_network->get_direction(aindex[0])[1] + forcey[0]*actin_network->get_direction(aindex[0])[0];
        force_par[1]=forcex[1]*actin_network->get_direction(aindex[1])[0] + forcey[1]*actin_network->get_direction(aindex[1])[1];
        force_perp[1]=-forcex[1]*actin_network->get_direction(aindex[1])[1] + forcey[1]*actin_network->get_direction(aindex[1])[0];

        torque[0]=cross(hx[0]-actin_network->get_position(aindex[0])[0],hy[0]-actin_network->get_position(aindex[0])[1],forcex[0],forcey[0]);
        torque[1]=cross(hx[1]-actin_network->get_position(aindex[1])[0],hy[1]-actin_network->get_position(aindex[1])[1],forcex[1],forcey[1]);
        actin_network->update_forces(aindex[0],force_par[0],force_perp[0],torque[0]);
        actin_network->update_forces(aindex[1],force_par[1],force_perp[1],torque[1]);
    }
    else
        return;

}

void motor::update_shape()
{
    if (state[0]==1 && state[1]==1) {
        hx[0]=actin_network->get_ends(aindex[0])[2]-pos_a_end[0]*actin_network->get_direction(aindex[0])[0];
        hy[0]=actin_network->get_ends(aindex[0])[3]-pos_a_end[0]*actin_network->get_direction(aindex[0])[1];
        //            pos_actin[0]=dis_points(hx[0],hy[0],actin_network->get_position(aindex[0])[0],actin_network->get_position(aindex[0])[1]);
        hx[1]=actin_network->get_ends(aindex[1])[2]-pos_a_end[1]*actin_network->get_direction(aindex[1])[0];
        hy[1]=actin_network->get_ends(aindex[1])[3]-pos_a_end[1]*actin_network->get_direction(aindex[1])[1];
        //            pos_actin[1]=dis_points(hx[1],hy[1],actin_network->get_position(aindex[1])[0],actin_network->get_position(aindex[1])[1]);
        mphi=atan2((hy[1]-hy[0]),(hx[1]-hx[0]));
    }
    else if(state[0]==1 && state[1]==0)
    {
        hx[0]=actin_network->get_ends(aindex[0])[2]-pos_a_end[0]*actin_network->get_direction(aindex[0])[0];
        hy[0]=actin_network->get_ends(aindex[0])[3]-pos_a_end[0]*actin_network->get_direction(aindex[0])[1];
        //            pos_actin[0]=dis_points(hx[0],hy[0],actin_network->get_position(aindex[0])[0],actin_network->get_position(aindex[0])[1]);
        mphi=atan2((hy[1]-hy[0]),(hx[1]-hx[0]));
    }
    else if(state[0]==0 && state[1]==1)
    {
        hx[1]=actin_network->get_ends(aindex[1])[2]-pos_a_end[1]*actin_network->get_direction(aindex[1])[0];
        hy[1]=actin_network->get_ends(aindex[1])[3]-pos_a_end[1]*actin_network->get_direction(aindex[1])[1];
        //            pos_actin[1]=dis_points(hx[1],hy[1],actin_network->get_position(aindex[1])[0],actin_network->get_position(aindex[1])[1]);
        mphi=atan2((hy[1]-hy[0]),(hx[1]-hx[0]));
    }
    else
    {
        return;
    }


}

inline void motor::move_end_detach(int hd, double speed, double pos)
{
    stretch=dis_points(hx[0],hy[0],hx[1],hy[1])-mld;
    std::cout<<"DEBUG: move_end_detach, before, hd = "<<hd<<" : color = "<< color<< "\tstretch = "<<stretch<<"\n";
    if (pos>=actin_network->get_alength(aindex[hd])) { //not sure why only "greater than or equal"
        if (event(kend,dt)==1) {
            std::cout<<"The new myosin position of head "<<hd<<" is OFF the actin filament AND detaching\n";
            state[hd]=0;
            aindex[hd]=-1;
            //                pos_actin[hd]=0;
            pos_a_end[hd]=0;
            hx[hd]=hx[pr(hd)]-pow(-1,hd)*mld*cos(mphi);
            hy[hd]=hy[pr(hd)]-pow(-1,hd)*mld*sin(mphi);
        }
        else {
            std::cout<<"The new myosin position of head "<<hd<<" is OFF the actin filament BUT not detaching\n";
            hx[hd]=actin_network->get_ends(aindex[hd])[2]-pos_a_end[hd]*actin_network->get_direction(aindex[hd])[0];
            hy[hd]=actin_network->get_ends(aindex[hd])[3]-pos_a_end[hd]*actin_network->get_direction(aindex[hd])[1];
            mphi=atan2((hy[1]-hy[0]),(hx[1]-hx[0]));
            //                pos_actin[hd]=dis_points(hx[hd],hy[hd],actin_network->get_position(aindex[hd])[0],actin_network->get_position(aindex[hd])[1]);
        }
    }
    else {
        std::cout<<"The new myosin position of head "<<hd<<" is "<<pos<<" away from the pointy end of the actin filament \n"; //still on the actin filament";
        pos_a_end[hd]=pos;
        hx[hd]=actin_network->get_ends(aindex[hd])[2]-pos_a_end[hd]*actin_network->get_direction(aindex[hd])[0];
        hy[hd]=actin_network->get_ends(aindex[hd])[3]-pos_a_end[hd]*actin_network->get_direction(aindex[hd])[1];
        mphi=atan2((hy[1]-hy[0]),(hx[1]-hx[0]));
        //            pos_actin[hd]=dis_points(hx[hd],hy[hd],actin_network->get_position(aindex[hd])[0],actin_network->get_position(aindex[hd])[1]);
    }
    stretch=dis_points(hx[0],hy[0],hx[1],hy[1])-mld;
    std::cout<<"DEBUG: move_end_detach, after, hd = "<<hd<<" : color = "<< color<< "\tstretch = "<<stretch<<"\n";

}

inline void motor::reflect(double x1, double x2, double y1, double y2)
{
    if (-fov[0]*0.5<x1<fov[0]*0.5 && -fov[0]*0.5<x2<fov[0]*0.5 && -fov[1]*0.5<y1<fov[1]*0.5 && -fov[1]*0.5<y2<fov[1]*0.5) {
        hx[0]=x1;
        hx[1]=x2;
        hy[0]=y1;
        hy[1]=y2;
    }
    else if (x1>=fov[0]*0.5 || x1<=-fov[0]*0.5)
    {
        hx[1]=x2;
        hy[0]=y1;
        hy[1]=y2;   
    }
    else if (x2>=fov[0]*0.5 || x2<=-fov[0]*0.5)
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

motor_ensemble::motor_ensemble(double mdensity, double fovx, double fovy, double mlen, actin_ensemble* network, double v0, double stiffness, double ron, double roff, double rend, double actin_len, double vis)
{
    fov[0]=fovx;
    fov[1]=fovy;
    mrho=mdensity;
    mld=mlen;
    nm=int(ceil(mrho*fov[0]*fov[1]));
    std::cout<<"\nDEBUG: Number of motors:"<<nm<<"\n";
    a_network=network;
    alpha=0.8;
    color = "0.5";//"green"; 
    for (int i=1; i<=nm; i++) {
        motorx=rng(-0.5*(fovx*alpha-mld),0.5*(fovx*alpha-mld));
        motory=rng(-0.5*(fovy*alpha-mld),0.5*(fovy*alpha-mld));
        mang=rng(0,2*pi);
        n_motors.push_back(motor(motorx,motory,mang,mld,a_network,0,0,-1,-1,fov[0],fov[1],v0,stiffness,ron,roff,rend,actin_len,vis,color));
    }
}

void motor_ensemble::motor_walk()
{

    for (int i=0; i<n_motors.size(); i++) {

        s[0]=n_motors[i].get_states()[0];
        s[1]=n_motors[i].get_states()[1];


        if (s[0]==0 && s[1]==0) {
            n_motors[i].attach(0);
            n_motors[i].attach(1);
            n_motors[i].brownian();
        }
        else if (s[0]==0 && s[1]==1) {
            n_motors[i].attach(0);
            n_motors[i].brownian();
            n_motors[i].step_onehead(1);
        }
        else if (s[0]==1 && s[1]==0) {
            n_motors[i].attach(1);
            n_motors[i].brownian();
            n_motors[i].step_onehead(0);
        }
        else {
            n_motors[i].step_twoheads();
        }

        n_motors[i].actin_update();
    }

}

void motor_ensemble::reshape()
{
    for (int i=0; i<n_motors.size(); i++) {
        n_motors[i].update_shape();
    }
}



void motor_ensemble::motor_write(std::ofstream& fout)
{
    for (int i=0; i<n_motors.size(); i++) {
        double stretch=dis_points(n_motors[i].get_heads()[0],n_motors[i].get_heads()[1],n_motors[i].get_heads()[2],n_motors[i].get_heads()[3])-mld;
        /*   if (stretch>3*0.25) {
             continue;
             }
             else{
             */   fout<<n_motors[i].get_heads()[0]<<"\t"<<n_motors[i].get_heads()[1]<<"\t"<<n_motors[i].get_heads()[2]-n_motors[i].get_heads()[0]<<"\t"<<n_motors[i].get_heads()[3]-n_motors[i].get_heads()[1]<<"\t"<<n_motors[i].get_color()<<"\n";
        //}
    } 
}

void motor_ensemble::motor_tension(std::ofstream& fout)
{
    for (int i=0; i<n_motors.size(); i++) {
        fout<<n_motors[i].tension()<<"\n";
    }
}

void motor_ensemble::add_motor(motor m)
{
    n_motors.push_back(m);
}
