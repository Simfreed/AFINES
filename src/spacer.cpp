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
    b_eng       = {0,0}; // filament / xlink bending energy
    
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
    b_eng       = {0,0}; // filament / xlink bending energy
    
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

void spacer::update_unattached()
{
    double new_rnd_x = rng_n(0,1), new_rnd_y = rng(0,1) new_rnd_th = rng_n(0,1);
    double w  = bd_rot_prefactor*(new_rnd_th +  prv_rnd_th);
    double vx = bd_prefactor*(new_rnd_x + prv_rnd_x);
    double vy = bd_prefactor*(new_rnd_y + prv_rnd_y);

    array<double, 2> pos = pos_bc(BC, actin_network->get_delrx(), dt, fov, {vx, vy}, {x + vx*dt, y + vy*dt});

    th += w*dt;
    x = pos[0];
    y = pos[1];

}

void spacer::update_head_pos_unattached()
{
    hx[0] = x - 0.5*mld*cos(th);
    hx[1] = x + 0.5*mld*cos(th);
    hy[0] = y - 0.5*mld*sin(th);
    hy[1] = y + 0.5*mld*sin(th)
}

//check for attachment of unbound heads given head index (0 for head 1, and 1 for head 2)
bool motor::attach(int hd)
{
    double not_off_prob = 0, onrate = kon;
    double mf_rand = rng(0,1.0);
    array<double, 2> intPoint;
    if (state[pr(hd)] == 1)
        onrate = kon2;
    
    set<pair<double, array<int, 2> > > dist_sorted = actin_network->get_binding_points(hx, hy, hd);

    if(!dist_sorted.empty()){
        

        for (set<pair<double, array<int, 2>>>::iterator it=dist_sorted.begin(); it!=dist_sorted.end(); ++it)
        {
            if (it->first > max_bind_dist) //since it's sorted, all the others will be farther than max_bind_dist too
                break;

            //head can't bind to the same filament link the other head is bound to
            else if(allowed_bind(hd, it->second)){
                
                intPoint = actin_network->get_filament((it->second).at(0))->get_link((it->second).at(1))->get_intpoint();
                not_off_prob += metropolis_prob(hd, it->second, intPoint, onrate);
                 
                if (mf_rand < not_off_prob) 
                {
                    //update state
                    state[hd] = 1;
                    f_index[hd] = (it->second).at(0);
                    this->set_l_index(hd, (it->second).at(1));
                    
                    //record displacement of head and orientation of link for future unbinding move
                    ldir_bind[hd] = actin_network->get_direction(f_index[hd], l_index[hd]);
                    bind_disp[hd] = rij_bc(BC, intPoint[0]-hx[hd], intPoint[1]-hy[hd], fov[0], fov[1], actin_network->get_delrx());

                    //update head position
                    hx[hd] = intPoint[0];
                    hy[hd] = intPoint[1];

                    //update relative head position
                    pos_a_end[hd]=dist_bc(BC, actin_network->get_end(f_index[hd], l_index[hd])[0] - hx[hd],
                                              actin_network->get_end(f_index[hd], l_index[hd])[1] - hy[hd], fov[0], fov[1], 
                                              actin_network->get_delrx());
                    
                    //(even if its at the barbed end upon binding, could have negative velocity, so always set this to false, until it steps)
                    at_barbed_end[hd] = false; 

                    return true;
                }
            }
        }
    }	
    return false;
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

/* ---------------------------------------------------------------------- */

void spacer::shake(int m)
{
  int nlist,list[2];
  double v[6]
  double invmass0,invmass1;

  // local atom IDs and constraint distances

  int i0 = atom->map(shake_atom[m][0]);
  int i1 = atom->map(shake_atom[m][1]);
  double bond1 = bond_distance[shake_type[m][0]];

  // r01 = distance vec between atoms, with PBC

  double r01[3];
  r01[0] = x[i0][0] - x[i1][0];
  r01[1] = x[i0][1] - x[i1][1];
  r01[2] = x[i0][2] - x[i1][2];
  domain->minimum_image(r01);

  // s01 = distance vec after unconstrained update, with PBC

  double s01[3];
  s01[0] = xshake[i0][0] - xshake[i1][0];
  s01[1] = xshake[i0][1] - xshake[i1][1];
  s01[2] = xshake[i0][2] - xshake[i1][2];
  domain->minimum_image(s01);

  // scalar distances between atoms

  double r01sq = r01[0]*r01[0] + r01[1]*r01[1] + r01[2]*r01[2];
  double s01sq = s01[0]*s01[0] + s01[1]*s01[1] + s01[2]*s01[2];

  // a,b,c = coeffs in quadratic equation for lamda

  if (rmass) {
    invmass0 = 1.0/rmass[i0];
    invmass1 = 1.0/rmass[i1];
  } else {
    invmass0 = 1.0/mass[type[i0]];
    invmass1 = 1.0/mass[type[i1]];
  }

  double a = (invmass0+invmass1)*(invmass0+invmass1) * r01sq;
  double b = 2.0 * (invmass0+invmass1) *
    (s01[0]*r01[0] + s01[1]*r01[1] + s01[2]*r01[2]);
  double c = s01sq - bond1*bond1;

  // error check

  double determ = b*b - 4.0*a*c;
  if (determ < 0.0) {
    error->warning(FLERR,"Shake determinant < 0.0",0);
    determ = 0.0;
  }

  // exact quadratic solution for lamda

  double lamda,lamda1,lamda2;
  lamda1 = (-b+sqrt(determ)) / (2.0*a);
  lamda2 = (-b-sqrt(determ)) / (2.0*a);

  if (fabs(lamda1) <= fabs(lamda2)) lamda = lamda1;
  else lamda = lamda2;

  // update forces if atom is owned by this processor

  lamda /= dtfsq;

  if (i0 < nlocal) {
    f[i0][0] += lamda*r01[0];
    f[i0][1] += lamda*r01[1];
    f[i0][2] += lamda*r01[2];
  }

  if (i1 < nlocal) {
    f[i1][0] -= lamda*r01[0];
    f[i1][1] -= lamda*r01[1];
    f[i1][2] -= lamda*r01[2];
  }

  if (evflag) {
    nlist = 0;
    if (i0 < nlocal) list[nlist++] = i0;
    if (i1 < nlocal) list[nlist++] = i1;

    v[0] = lamda*r01[0]*r01[0];
    v[1] = lamda*r01[1]*r01[1];
    v[2] = lamda*r01[2]*r01[2];
    v[3] = lamda*r01[0]*r01[1];
    v[4] = lamda*r01[0]*r01[2];
    v[5] = lamda*r01[1]*r01[2];

    v_tally(nlist,list,2.0,v);
  }
}
