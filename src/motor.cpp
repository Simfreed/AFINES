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
#include "filament_ensemble.h"
//#include "actin.h"

//motor class
template <class filament_ensemble_type>
motor<filament_ensemble_type>::motor( array<double, 3> pos, double mlen, filament_ensemble_type * network, 
        array<int, 2> mystate, array<int, 2> myfindex, array<int, 2> mylindex,
        array<double, 2> myfov, double delta_t, double temp,
        double v0, double stiffness, double ron, double roff, double
        rend, double actin_len, double vis, string bc) {
    
    vs          = v0;//rng_n(v0,0.4);//rng(v0-0.3,v0+0.3);
    dm          = 0.25;//actin_len/10; //max binding distance
    mk          = stiffness;//rng(10,100); 
    fmax        = mk*dm*2;//rng(1,20);
    mld         = mlen;
    kon         = ron;
    koff        = roff;
    kend        = rend;
    mphi        = pos[2];
    dt          = delta_t;
    temperature = temp;
    state       = mystate;
    f_index     = myfindex; //filament index for each head
    l_index     = mylindex; //link index for each head
    fov         = myfov;
    BC          = bc; 
    
    actin_network = network;
    
    hx[0]=pos[0]-0.5*mld*cos(mphi);
    hy[0]=pos[1]-0.5*mld*sin(mphi);
    hx[1]=hx[0]+mld*cos(mphi);
    hy[1]=hy[0]+mld*sin(mphi);
    
    mobility=log(10)/(4*pi*vis*mld);
    
    
    // pos_a_end = distance from pointy end -- by default 0
    //      i.e., if l_index[hd] = j, then pos_a_end[hd] is the distance to the "j+1"th actin
    pos_a_end = {0, 0};

    if (state[0]){
        pos_a_end[0] = dis_points(hx[0],hy[0],
                actin_network->get_end(f_index[0], l_index[0])[0],
                actin_network->get_end(f_index[0], l_index[0])[1]);
    }
    if (state[1]){
        pos_a_end[1] = dis_points(hx[1],hy[1],
                actin_network->get_end(f_index[1], l_index[1])[0],
                actin_network->get_end(f_index[1], l_index[1])[1]);
    }
    

}

template <class filament_ensemble_type>
motor<filament_ensemble_type>::~motor(){ 
};

//return motor state with a given head number
template <class filament_ensemble_type>
array<int, 2> motor<filament_ensemble_type>::get_states() 
{
    return state;
}

template <class filament_ensemble_type>
array<double, 2> motor<filament_ensemble_type>::get_hx()
{
    return hx;
}

template <class filament_ensemble_type>
array<double, 2> motor<filament_ensemble_type>::get_hy()
{
    return hy;
}

template <class filament_ensemble_type>
string motor<filament_ensemble_type>::get_BC()
{
    return BC;
}

template <class filament_ensemble_type>
double motor<filament_ensemble_type>::tension()
{
    double lf=dis_points(hx[0],hy[0],hx[1],hy[1]);
    return mk*(lf-mld)/fmax;
}

//check for attachment of unbound heads given head index (0 for head 1, and 1 for head 2)
template <class filament_ensemble_type>
void motor<filament_ensemble_type>::attach(int hd)
{
    map<array<int, 2>, double> dist = actin_network->get_dist(hx[hd],hy[hd]);
    double onrate;
    array<double, 2> intpoint;

    if(!dist.empty()){
        for (map<array<int, 2>, double>::iterator it=dist.begin(); it!=dist.end(); ++it)
        { 
            if (it->second <= dm && f_index[pr(hd)]!=(it->first).at(0) && l_index[pr(hd)] != (it->first).at(1)) {
                
                onrate=kon*exp(-((it->second)*(it->second))/(dm*dm));
                
                if (event(onrate,dt)==1) {
                    
                    //update state
                    state[hd] = 1;
                    f_index[hd] = (it->first).at(0);
                    l_index[hd] = (it->first).at(1);

                    //update head position
                    intpoint = actin_network->get_intpoints(f_index[hd], l_index[hd], hx[hd],hy[hd]);
                    hx[hd] = intpoint[0];
                    hy[hd] = intpoint[1];

                    if (state[pr(hd)]==0) {
                        hx[pr(hd)] = hx[hd]+pow(-1,hd)*mld*cos(mphi);
                        hy[pr(hd)] = hy[hd]+pow(-1,hd)*mld*sin(mphi);
                    }
                    else {
                        mphi=atan2((hy[1]-hy[0]),(hx[1]-hx[0]));
                    }

                    pos_a_end[hd]=dis_points(hx[hd],hy[hd],
                            actin_network->get_end(f_index[hd], l_index[hd])[0],
                            actin_network->get_end(f_index[hd], l_index[hd])[1]);
                    break;
                }
            }
        }
    }	
} 

//perform brownian motion and shear if head unattached
template <class filament_ensemble_type>
void motor<filament_ensemble_type>::brownian(double t, double gamma)
{
    if (state[0]==0 && state[1]==0) {

        xm[0]=hx[0]+sqrt(dt*mobility*2*temperature)*rng_n(0,1) - mk*dt*(hx[0]-hx[1]+mld*cos(mphi)) + gamma*dt*hy[0];
        xm[1]=hx[1]+sqrt(dt*mobility*2*temperature)*rng_n(0,1) + mk*dt*(hx[0]-hx[1]+mld*cos(mphi)) + gamma*dt*hy[1];
        ym[0]=hy[0]+sqrt(dt*mobility*2*temperature)*rng_n(0,1) - mk*dt*(hy[0]-hy[1]+mld*sin(mphi));
        ym[1]=hy[1]+sqrt(dt*mobility*2*temperature)*rng_n(0,1) + mk*dt*(hy[0]-hy[1]+mld*sin(mphi));
        if (BC == "REFLECTIVE") reflect(t, gamma, xm[0],xm[1],ym[0],ym[1]);
        if (BC == "PERIODIC")  periodic(t, gamma, xm[0],xm[1],ym[0],ym[1]);
        mphi=atan2((hy[1]-hy[0]),(hx[1]-hx[0]));
    }
    else if (state[0]==0 || state[1]==0) {
        int hd=state[0];//cleverly equivalent to: 
                        //  int hd = (hd for which state[hd] = 0)

        xm[hd]=hx[hd]+sqrt(dt*mobility*2*temperature)*rng_n(0,1) - mk*dt*(hx[hd]-hx[pr(hd)]+pow(-1,hd)*mld*cos(mphi)) + gamma*dt*hy[hd];
        ym[hd]=hy[hd]+sqrt(dt*mobility*2*temperature)*rng_n(0,1) - mk*dt*(hy[hd]-hy[pr(hd)]+pow(-1,hd)*mld*sin(mphi));
        xm[pr(hd)]=hx[pr(hd)];
        ym[pr(hd)]=hy[pr(hd)];
        if (BC == "REFLECTIVE") reflect(t, gamma, xm[0],xm[1],ym[0],ym[1]);
        if (BC == "PERIODIC")  periodic(t, gamma, xm[0],xm[1],ym[0],ym[1]);
        mphi=atan2((hy[1]-hy[0]),(hx[1]-hx[0]));
    }
    else {
        return;
    }

}

//stepping and detachment kinetics of a single bound head 
template <class filament_ensemble_type>
void motor<filament_ensemble_type>::step_onehead(int hd)
{

    if (event(koff,dt)==1) {
        this->detach_head(hd);
    }
    else
    {
        move_end_detach(hd, pos_a_end[hd]+dt*vs);
    }
}

//stepping and detachment kinetics of a motor with both heads attached
template <class filament_ensemble_type>
void motor<filament_ensemble_type>::step_twoheads()
{
    array<double, 2> vm, fm, offrate;
        
    //fm = vec(Fm).(-vec(u)) 
    fm[0]=mk*((hx[0]-hx[1]+mld*cos(mphi))*actin_network->get_direction(f_index[0],l_index[0])[0] +
              (hy[0]-hy[1]+mld*sin(mphi))*actin_network->get_direction(f_index[0],l_index[0])[1]); 
    fm[1]=mk*(-(hx[0]-hx[1]+mld*cos(mphi))*actin_network->get_direction(f_index[1],l_index[1])[0] -
            (hy[0]-hy[1]+mld*sin(mphi))*actin_network->get_direction(f_index[1],l_index[1])[1]);

    
    vm[0]=velocity(vs,fm[0],fmax);
    vm[1]=velocity(vs,fm[1],fmax);
    
    /*impose force-dependent bell's law on detachment rates*/
    offrate[0]=koff*exp(fabs(fm[0])/fmax);
    offrate[1]=koff*exp(fabs(fm[1])/fmax);

    if (event(offrate[0],dt)==1) {
        this->detach_head(0);
        move_end_detach(1, pos_a_end[1]+dt*vs);
    }

    else if (event(offrate[1],dt)==1) {
        this->detach_head(1);
        move_end_detach(0,pos_a_end[0]+dt*vs);
    }

    else {
        move_end_detach(0,pos_a_end[0]+dt*vm[0]);
        move_end_detach(1,pos_a_end[1]+dt*vm[1]); 
    }
}

// Using the lever rule to propagate force as outlined in Nedelec F 2002
template <class filament_ensemble_type>
void motor<filament_ensemble_type>::actin_update()
{
    if (state[0]==1 && state[1]==1) {

        array<double, 2> fx, fy, pos_ratio;
        
        fx[0]= -mk*(hx[0]-hx[1]+mld*cos(mphi));
        fx[1]= -fx[0];
        fy[0]= -mk*(hy[0]-hy[1]+mld*sin(mphi));
        fy[1]= -fy[0];
        pos_ratio[0] = pos_a_end[0]/actin_network->get_llength(f_index[0], l_index[0]);
        pos_ratio[1] = pos_a_end[1]/actin_network->get_llength(f_index[1], l_index[1]);

        actin_network->update_forces(f_index[0], l_index[0] + 1, fx[0] *    pos_ratio[0] ,  fy[0] *    pos_ratio[0]);
        actin_network->update_forces(f_index[0], l_index[0]    , fx[0] * (1-pos_ratio[0]), fy[0] * (1-pos_ratio[0]));
        actin_network->update_forces(f_index[1], l_index[1] + 1, fx[1] *    pos_ratio[1] , fy[1] *    pos_ratio[1]);
        actin_network->update_forces(f_index[1], l_index[1]    , fx[1] * (1-pos_ratio[1]), fy[1] * (1-pos_ratio[1]));
        
    }
    else
        return;

}

template <class filament_ensemble_type>
void motor<filament_ensemble_type>::update_shape()
{
    if (state[0]==1 && state[1]==1) {
        hx[0]=actin_network->get_end(f_index[0],l_index[0])[0]-pos_a_end[0]*actin_network->get_direction(f_index[0],l_index[0])[0];
        hy[0]=actin_network->get_end(f_index[0],l_index[0])[1]-pos_a_end[0]*actin_network->get_direction(f_index[0],l_index[0])[1];
        hx[1]=actin_network->get_end(f_index[1],l_index[1])[0]-pos_a_end[1]*actin_network->get_direction(f_index[1],l_index[1])[0];
        hy[1]=actin_network->get_end(f_index[1],l_index[1])[1]-pos_a_end[1]*actin_network->get_direction(f_index[1],l_index[1])[1];
        mphi=atan2((hy[1]-hy[0]),(hx[1]-hx[0]));
    }
    else if(state[0]==1 && state[1]==0)
    {
        hx[0]=actin_network->get_end(f_index[0],l_index[0])[0]-pos_a_end[0]*actin_network->get_direction(f_index[0],l_index[0])[0];
        hy[0]=actin_network->get_end(f_index[0],l_index[0])[1]-pos_a_end[0]*actin_network->get_direction(f_index[0],l_index[0])[1];
        mphi=atan2((hy[1]-hy[0]),(hx[1]-hx[0]));
    }
    else if(state[0]==0 && state[1]==1)
    {
        hx[1]=actin_network->get_end(f_index[1],l_index[1])[0]-pos_a_end[1]*actin_network->get_direction(f_index[1],l_index[1])[0];
        hy[1]=actin_network->get_end(f_index[1],l_index[1])[1]-pos_a_end[1]*actin_network->get_direction(f_index[1],l_index[1])[1];
        mphi=atan2((hy[1]-hy[0]),(hx[1]-hx[0]));
    }
    else
    {
        return;
    }

}

template <class filament_ensemble_type>
void motor<filament_ensemble_type>::detach_head(int hd)
{
    state[hd]=0;
    f_index[hd]=-1;
    l_index[hd]=-1;

    pos_a_end[hd]=0;
    hx[hd]=hx[pr(hd)]-pow(-1,hd)*mld*cos(mphi);
    hy[hd]=hy[pr(hd)]-pow(-1,hd)*mld*sin(mphi);

}

template <class filament_ensemble_type>
void motor<filament_ensemble_type>::move_end_detach(int hd, double pos)
{
    double link_length = actin_network->get_llength(f_index[hd],l_index[hd]);
    if (pos >= link_length) { // "passed" the link
        
        if (l_index[hd] == 0){ // the barbed end of the filament
            
            if (event(kend,dt)==1) {
 
                this->detach_head(hd);

            }
            else {

                hx[hd]=actin_network->get_end(f_index[hd],l_index[hd])[0]-pos_a_end[hd]*actin_network->get_direction(f_index[hd],l_index[hd])[0];
                hy[hd]=actin_network->get_end(f_index[hd],l_index[hd])[1]-pos_a_end[hd]*actin_network->get_direction(f_index[hd],l_index[hd])[1];
                mphi=atan2((hy[1]-hy[0]),(hx[1]-hx[0]));
            
            }
        }
        else{ 
            /*Move the motor to the next link on the filament
             *At the projected new position along that filament
             */
            
            l_index[hd] = l_index[hd] - 1;
            pos_a_end[hd] = pos - link_length;
            hx[hd]=actin_network->get_end(f_index[hd],l_index[hd])[0]-pos_a_end[hd]*actin_network->get_direction(f_index[hd],l_index[hd])[0];
            hy[hd]=actin_network->get_end(f_index[hd],l_index[hd])[1]-pos_a_end[hd]*actin_network->get_direction(f_index[hd],l_index[hd])[1];
            mphi=atan2((hy[1]-hy[0]),(hx[1]-hx[0]));
    
        }
    }
    else {
        pos_a_end[hd]=pos;
        hx[hd]=actin_network->get_end(f_index[hd],l_index[hd])[0]-pos_a_end[hd]*actin_network->get_direction(f_index[hd],l_index[hd])[0];
        hy[hd]=actin_network->get_end(f_index[hd],l_index[hd])[1]-pos_a_end[hd]*actin_network->get_direction(f_index[hd],l_index[hd])[1];
        mphi=atan2((hy[1]-hy[0]),(hx[1]-hx[0]));
    }
    

}

template <class filament_ensemble_type>
inline void motor<filament_ensemble_type>::reflect(double t, double gamma, double x1, double x2, double y1, double y2)
{
    //Calculate the sheared simulation bounds (at this height)
    double xleft = 0, xright = 0, yleft, yright;
    xleft  =  max(-fov[0] * 0.5 + gamma * y1 * t, -fov[0] * 0.5 + gamma * y2 * t);
    xright =  min( fov[0] * 0.5 + gamma * y1 * t,  fov[0] * 0.5 + gamma * y2 * t);
    yleft  = -fov[1]*0.5;
    yright =  fov[1]*0.5;
    
    if (    xleft < x1 && x1 < xright &&  xleft < x2 && x2 < xright
        &&  yleft < y1 && y1 < yright &&  yleft < y2 && y2 < yright) {
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

//TODO: Implement harmonic boundary conditions
template <class filament_ensemble_type>
inline void motor<filament_ensemble_type>::periodic(double t, double gamma, double x1, double x2, double y1, double y2)
{
    //Calculate the sheared simulation bounds (at this height)
    double xleft, xright, yleft, yright;
    xleft  =  max(-fov[0] * 0.5 + gamma * y1 * t, -fov[0] * 0.5 + gamma * y2 * t);
    xright =  min( fov[0] * 0.5 + gamma * y1 * t,  fov[0] * 0.5 + gamma * y2 * t);
    yleft  = -fov[1]*0.5;
    yright =  fov[1]*0.5;
    
    if (    xleft < x1 && x1 < xright &&  xleft < x2 && x2 < xright
        &&  yleft < y1 && y1 < yright &&  yleft < y2 && y2 < yright) {
        hx[0]=x1;
        hx[1]=x2;
        hy[0]=y1;
        hy[1]=y2;
    }
    else if (x1>=xright || x1<=xleft)
    {
        if( x1 >= xright) hx[0]=x1 - fov[0];
        if( x1 <= xleft ) hx[0]=x1 + fov[0];
        hx[1]=x2;
        hy[0]=y1; 
        hy[1]=y2;   
    }
    else if (x2>=xright || x2<=xleft)
    {
        hx[0]=x1;
        if( x2 >= xright) hx[1]=x2 - fov[0];
        if( x2 <= xleft ) hx[1]=x2 + fov[0];
        hy[0]=y1;
        hy[1]=y2;
    }
    else if(y1>=yright || y1<=yleft)
    {
        hx[0]=x1;
        hx[1]=x2;
        if( y1 >= yright) hy[0]=y1 - fov[1];
        if( y1 <= yleft ) hy[0]=y1 + fov[1];
        hy[1]=y2;
    }
    else{
        hx[0]=x1;
        hx[1]=x2;
        hy[0]=y1;
        if( y2 >= yright) hy[1]=y2 - fov[1];
        if( y2 <= yleft ) hy[1]=y2 + fov[1];
    }
}

template <class filament_ensemble_type>
array<int, 2> motor<filament_ensemble_type>::get_f_index(){
    return f_index;
}

template <class filament_ensemble_type>
array<int, 2> motor<filament_ensemble_type>::get_l_index(){
    return l_index;
}

template <class filament_ensemble_type>
array<double, 2> motor<filament_ensemble_type>::get_pos_a_end(){
    return pos_a_end;
}

template class motor<ATfilament_ensemble>;
template class motor<baoab_filament_ensemble>;
template class motor<lammps_filament_ensemble>;
template class motor<langevin_leapfrog_filament_ensemble>;
