/*
 * spacer.cpp
 *  
 *
 *  Created by Simon Freedman
 *  Copyright 2013 University of Chicago. All rights reserved.
 *
 */

#include "spacer.h"
#include "filament_ensemble.h"
//spacer class

spacer::spacer( array<double, 3> pos, 
        double mlen, filament_ensemble * network, 
        array<int, 2> mystate, 
        array<int, 2> myfindex, 
        array<int, 2> mylindex,
        array<double, 2> myfov, 
        double delta_t, 
        double temp,
        double v0, 
        double stiffness, 
        double max_ext_ratio, 
        double ron, double roff, double rend, 
        double fstall, double rcut,
        double vis, string bc) {
    
    vs          = v0;
    mk          = stiffness;//rng(10,100); 
    
    stall_force   = fstall;
    temperature   = temp;

    max_bind_dist = rcut;

    mld         = mlen;
    dt          = delta_t;
    kon         = ron*dt;
    koff        = roff*dt;
    kend        = rend*dt;
    kon2        = ron*dt;
    koff2       = roff*dt;
    kend2       = rend*dt;
    mphi        = pos[2];
    state       = mystate;
    f_index     = myfindex; //filament index for each head
    l_index     = mylindex; //link index for each head
    fov         = myfov;
    BC          = bc; 
    actin_network = network;
    damp        =(6*pi*vis*mld);
    bd_prefactor= sqrt(temperature/(2*damp*dt)); 
    
    /****for FENE motors******/
    max_ext     = max_ext_ratio*mlen;
    eps_ext     = 0.01*max_ext;
    /*************************/
    
    shear       = 0;
    tension     = 0;
    force       = {0,0}; // tensile force on the spring  
    
    b_force[0]  = {0,0}; //b_force[0] = bending force on head 0 due to h0-h1-link angle in cartesian coords
    b_force[1]  = {0,0}; //b_force[1] = bending force on head 1 due to h1-h0-link angle in cartesian coords

    kinetic_energy = 0; //assume m = 1
    
    array<double, 2> posH0 = boundary_check(0, pos[0]-0.5*mld*cos(mphi), pos[1]-0.5*mld*sin(mphi)); 
    array<double, 2> posH1 = boundary_check(1, pos[0]+0.5*mld*cos(mphi), pos[1]+0.5*mld*sin(mphi)); 
    hx[0] = posH0[0];
    hy[0] = posH0[1];
    hx[1] = posH1[0];
    hy[1] = posH1[1];
    
    disp = rij_bc(BC, hx[1]-hx[0], hy[1]-hy[0], fov[0], fov[1], actin_network->get_delrx()); 
    
    pos_a_end = {0, 0}; // pos_a_end = distance from pointy end -- by default 0
                        // i.e., if l_index[hd] = j, then pos_a_end[hd] is the distance to the "j+1"th actin
    
    ldir_bind[0] = {0,0};
    ldir_bind[1] = {0,0};

    bind_disp[0] = {0,0};
    bind_disp[1]=  {0,0};

    at_barbed_end = {false, false};
    if (state[0] == 1){
        pos_a_end[0] = dist_bc(BC, actin_network->get_end(f_index[0], l_index[0])[0] - hx[0],
                actin_network->get_end(f_index[0], l_index[0])[1] - hy[0], fov[0], fov[1], 0);
        ldir_bind[0] = actin_network->get_direction(f_index[0], l_index[0]);

    }
    if (state[1] == 1){
        pos_a_end[1] = dist_bc(BC, actin_network->get_end(f_index[1], l_index[1])[0] - hx[1],
                                   actin_network->get_end(f_index[1], l_index[1])[1] - hy[1], fov[0], fov[1], 0);
        ldir_bind[1] = actin_network->get_direction(f_index[1], l_index[1]);
    }
    
    prv_rnd_x = {0,0};
    prv_rnd_y = {0,0};

}


spacer::spacer( array<double, 4> pos, 
        double mlen, filament_ensemble * network, 
        array<int, 2> mystate, 
        array<int, 2> myfindex, 
        array<int, 2> mylindex,
        array<double, 2> myfov, 
        double delta_t, 
        double temp,
        double v0, 
        double stiffness, 
        double max_ext_ratio, 
        double ron, double roff, double rend, 
        double fstall, double rcut,
        double vis, string bc) {
    
    vs          = v0;
    mk          = stiffness;
    
    stall_force = fstall;
    temperature = temp;

    max_bind_dist = rcut;
    
    mld         = mlen;
    dt          = delta_t;
    kon         = ron*dt;
    koff        = roff*dt;
    kend        = rend*dt;
    kon2        = ron*dt;
    koff2       = roff*dt;
    kend2       = rend*dt;
    mphi        = pos[2];
    state       = mystate;
    f_index     = myfindex; //filament index for each head
    l_index     = mylindex; //link index for each head
    fov         = myfov;
    BC          = bc; 
    actin_network = network;
    damp        =(6*pi*vis*mld);
    bd_prefactor= sqrt(temperature/(2*damp*dt)); 
    
    /********for FENE springs*********/ 
    max_ext     = max_ext_ratio*mlen;
    eps_ext     = 0.01*max_ext;
    /********************************/

    shear       = 0;
    tension     = 0;
    force       = {0,0}; // force on the spring  
    
    b_force[0]  = {0,0}; //b_force[0] = bending force on head 0 due to h0-h1-link angle in cartesian coords
    b_force[1]  = {0,0}; //b_force[1] = bending force on head 1 due to h1-h0-link angle in cartesian coords
    
    kinetic_energy = 0;
    pos_a_end = {0, 0}; // pos_a_end = distance from pointy end -- by default 0
                        // i.e., if l_index[hd] = j, then pos_a_end[hd] is the distance to the "j+1"th actin

    array<double, 2> posH0 = boundary_check(0, pos[0], pos[1]); 
    array<double, 2> posH1 = boundary_check(1, pos[0]+pos[2], pos[1]+pos[3]); 
    hx[0] = posH0[0];
    hy[0] = posH0[1];
    hx[1] = posH1[0];
    hy[1] = posH1[1];
   
    //force can be non-zero and angle is determined from disp vector
    this->update_angle();
    this->update_force();
    
    ldir_bind[0] = {0,0};
    ldir_bind[1] = {0,0};
    bind_disp[0] = {0,0};
    bind_disp[1] = {0,0};

    at_barbed_end = {false, false};

    if (state[0] == 1){
        pos_a_end[0] = dist_bc(BC, actin_network->get_end(f_index[0], l_index[0])[0] - hx[0],
                                   actin_network->get_end(f_index[0], l_index[0])[1] - hy[0], fov[0], fov[1], 0);
        ldir_bind[0] = actin_network->get_direction(f_index[0], l_index[0]);

    }
    if (state[1] == 1){
        pos_a_end[1] = dist_bc(BC, actin_network->get_end(f_index[1], l_index[1])[0] - hx[1],
                                   actin_network->get_end(f_index[1], l_index[1])[1] - hy[1], fov[0], fov[1], 0);
        ldir_bind[1] = actin_network->get_direction(f_index[1], l_index[1]);
    }

    prv_rnd_x = {0,0};
    prv_rnd_y = {0,0};

}

 spacer::~spacer(){};

void spacer::set_bending(double force_constant, double ang){
    kb  = force_constant;
    th0 = ang;
}

void spacer::update_force()
{ 
    //cout<<"\nDEBUG: using spacer update_force()";
    if (state[0] == 1 && state[1] == 1){
        update_bending(0);
        update_bending(1);
    }
    else
    {
        b_force[0] = {0,0};
        b_force[1] = {0,0};
    }

    update_angle(); // need to recalculate this, since heads might've moved

    tension = (mk)*(hypot(disp[0], disp[1]) - mld);
    force = {tension*cos(mphi), tension*sin(mphi)};
    
}
  //Measure distance to FARTHER END of actin filament that the spacer is bound to
  //
int spacer::get_further_end(int hd, int findex, int lindex)
{
    return (pos_a_end[hd] > 0.5*actin_network->get_llength(findex, lindex));
}

array<double, 2> spacer::disp_from_actin(int hd, int findex, int aindex)
{
  return rij_bc(BC,
          actin_network->get_filament(findex)->get_actin(aindex)->get_xcm() - hx[hd],
          actin_network->get_filament(findex)->get_actin(aindex)->get_ycm() - hy[hd],
          fov[0], fov[1], actin_network->get_delrx());
}

void spacer::update_bending(int hd)
{
  array<double, 2> delr1, delr2;
  double f1[2], f3[2];
  double rsq1,rsq2,r1,r2,c,s,dth,a,a11,a12,a22;
  
  int actin_further_end = get_further_end(hd, f_index[hd], l_index[hd]);

  // 1st bond
  delr1 = disp_from_actin(hd, f_index[hd], l_index[hd] + actin_further_end); 
  rsq1  = delr1[0]*delr1[0] + delr1[1]*delr1[1];
  r1    = sqrt(rsq1);

  // 2nd bond
  delr2 = {pow(-1, hd)*disp[0], pow(-1, hd)*disp[1]};
  rsq2  = delr2[0]*delr2[0] + delr2[1]*delr2[1];
  r2    = sqrt(rsq2);

  // angle (cos and sin)
  c = (delr1[0]*delr2[0] + delr1[1]*delr2[1]) / (r1*r2);
    
  if (c > 1.0) c = 1.0;
  if (c < -1.0) c = -1.0;

  s = sqrt(1.0 - c*c);
  if (s < maxSmallAngle) s = maxSmallAngle;

  dth = acos(c) - th0;

  // force
  a   = -kb * dth / s; 
  a11 = a*c / rsq1;
  a12 = -a / (r1*r2);
  a22 = a*c / rsq2;

  f1[0] = a11*delr1[0] + a12*delr2[0];
  f1[1] = a11*delr1[1] + a12*delr2[1];
  f3[0] = a22*delr2[0] + a12*delr1[0];
  f3[1] = a22*delr2[1] + a12*delr1[1];

  // apply force to each of 3 atoms
  
  actin_network->update_forces(f_index[hd], l_index[hd] + actin_further_end, f1[0], f1[1]);
  b_force[hd][0] += (-f1[0] - f3[0]);
  b_force[hd][1] += (-f1[1] - f3[1]);
  b_force[pr(hd)][0] += f3[0];
  b_force[pr(hd)][1] += f3[1];
  
  b_eng[hd] = kb*dth*dth/(r1+r2);

}
double spacer::get_kb(){
    return kb;
}

double spacer::get_th0(){
    return th0;
}

void spacer::identify(){
    cout<<"I am a spacer";
}

void spacer::brownian_relax(int hd)
{
    //cout<<"\nDEBUG: using spacer brownian_relax()";
    
    double new_rnd_x= rng_n(0,1), new_rnd_y = rng_n(0,1);
    
    double vx =  (pow(-1,hd)*force[0] + b_force[hd][0]) / damp + bd_prefactor*(new_rnd_x + prv_rnd_x[hd]);
    double vy =  (pow(-1,hd)*force[1] + b_force[hd][1]) / damp + bd_prefactor*(new_rnd_y + prv_rnd_y[hd]);
    kinetic_energy = vx*vx + vy*vy;    
    array<double, 2> pos = boundary_check(hd, hx[hd] + vx*dt, hy[hd] + vy*dt);
    hx[hd] = pos[0];
    hy[hd] = pos[1];

    prv_rnd_x[hd] = new_rnd_x;
    prv_rnd_y[hd] = new_rnd_y;

}

void spacer::actin_update()
{
    if (state[0]==1) this->actin_update_hd(0, { force[0] + b_force[0][0],  force[1] + b_force[0][1]});
    if (state[1]==1) this->actin_update_hd(1, {-force[0] + b_force[1][0], -force[1] + b_force[1][1]});
    
    //reset bending force
    b_force[0] = {0,0};
    b_force[1] = {0,0};
}

array<array<double, 2>,2> spacer::get_b_force()
{
    return b_force;
}


//metropolis algorithm with rate constant
double spacer::metropolis_prob(int hd, array<int, 2> fl_idx, array<double, 2> newpos, double maxprob)
{
//    cout<<"\nDEBUG: using spacer metropolis_prob";
    double prob = maxprob;
    double stretch  = dist_bc(BC, newpos[0] - hx[pr(hd)], newpos[1] - hy[pr(hd)], fov[0], fov[1], actin_network->get_delrx()) - mld; 
    
    array<double, 2> delr1, delr2;
    double r1, r2, c, dth, bend_eng = 0;
    
    if (state[hd] == 0 && state[pr(hd)] == 1) { //it's trying to attach
        
//        delr1 = disp_from_actin(hd, it->second.at(0), it->second.at(1) + get_further_end(hd, it->second.at(0), it->second.at(1))); 
        delr1 = disp_from_actin(hd, fl_idx[0], fl_idx[1] + get_further_end(hd, fl_idx[0], fl_idx[1])); 
        r1  = sqrt(delr1[0]*delr1[0] + delr1[1]*delr1[1]);

        // 2nd bond
        delr2 = {pow(-1, hd)*disp[0], pow(-1, hd)*disp[1]};
        r2  = sqrt(delr2[0]*delr2[0] + delr2[1]*delr2[1]);

        // cos
        c = (delr1[0]*delr2[0] + delr1[1]*delr2[1]) / (r1*r2);

        if (c > 1.0) c = 1.0;
        if (c < -1.0) c = -1.0;

        dth = acos(c) - th0;
        bend_eng = kb*dth*dth/(r1+r2);
        
    }
    
    double delE = 0.5*mk*stretch*stretch + bend_eng - this->get_stretching_energy() - b_eng[hd];

    if( delE > 0 )
        prob *= exp(-delE/temperature);
    
    return prob;
}

bool spacer::allowed_bind(int hd, array<int, 2> fl_idx)
{
//    cout<<"\nDEBUG: using spacer allowed bind";
    return (fl_idx[0] != f_index[pr(hd)]);
}
