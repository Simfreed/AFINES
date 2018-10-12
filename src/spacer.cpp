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

    max_bind_dist    = rcut;
    max_bind_dist_sq = rcut*rcut;

    mld         = mlen;
    dt          = delta_t;
    kon         = ron*dt;
    koff        = roff*dt;
    kend        = rend*dt;
    kon2        = ron*dt;
    koff2       = roff*dt;
    kend2       = rend*dt;
    state       = mystate;
    f_index     = myfindex; //filament index for each head
    
    filament_network = network;
    init_l_index(0, mylindex[0]);
    init_l_index(1, mylindex[1]);
    
    fov         = myfov;
    BC          = bc; 
    damp        =(6*pi*vis*mld);
    bd_prefactor= sqrt(temperature/(2*damp*dt)); 
    rot_bd_prefactor= sqrt(temperature/(2*damp*dt)); 

    /****for FENE motors******/
    max_ext     = max_ext_ratio*mlen;
    eps_ext     = 0.01*max_ext;
    /*************************/

    shear       = 0;
    tension     = 0;
    force       = {{0,0}}; // tensile force on the spring  
    b_eng       = {{0,0}}; // filament / xlink bending energy
    
    b_force[0]  = {{0,0}}; //b_force[0] = bending force on head 0 due to h0-h1-spring angle in cartesian coords
    b_force[1]  = {{0,0}}; //b_force[1] = bending force on head 1 due to h1-h0-spring angle in cartesian coords

    kinetic_energy = 0; //assume m = 1

    cm = {{pos[0], pos[1]}};
    th = pos[2];
    array<double, 2> posH0 = boundary_check(0, pos[0]-0.5*mld*cos(pos[2]), pos[1]-0.5*mld*sin(pos[2])); 
    array<double, 2> posH1 = boundary_check(1, pos[0]+0.5*mld*cos(pos[2]), pos[1]+0.5*mld*sin(pos[2])); 
    hx[0] = posH0[0];
    hy[0] = posH0[1];
    hx[1] = posH1[0];
    hy[1] = posH1[1];
    
    disp = rij_bc(BC, hx[1]-hx[0], hy[1]-hy[0], fov[0], fov[1], filament_network->get_delrx()); 
    
    pos_rat = {{0, 0}}; // pos_rat = distance from pointy end / segment length -- by default 0
                        // i.e., if l_index[hd] = j, then pos_rat[hd]*|disp[j]| is the distance to the "j+1"th bead
    
    ldir_bind[0] = {{0,0}};
    ldir_bind[1] = {{0,0}};

    bind_disp[0] = {{0,0}};
    bind_disp[1]=  {{0,0}};
    
    bind_rot = {{0,0}};

    at_barbed_end = {{false, false}};

    if (state[0] == 1){
        update_pos_rat(0);
        ldir_bind[0] = filament_network->get_direction(f_index[0], l_index[0]);

    }
    if (state[1] == 1){
        update_pos_rat(1);
        ldir_bind[1] = filament_network->get_direction(f_index[1], l_index[1]);
    }
    
    prv_rnd_x = 0;
    prv_rnd_y = 0;
    prv_rnd_th = 0;
    
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
    max_bind_dist_sq = rcut*rcut;
    
    mld         = mlen;
    dt          = delta_t;
    kon         = ron*dt;
    koff        = roff*dt;
    kend        = rend*dt;
    kon2        = ron*dt;
    koff2       = roff*dt;
    kend2       = rend*dt;

    state       = mystate;
    f_index     = myfindex; //filament index for each head
    filament_network = network;
    init_l_index(0, mylindex[0]);
    init_l_index(1, mylindex[1]);
    fov         = myfov;
    BC          = bc; 
    damp        =(6*pi*vis*mld);
    bd_prefactor= sqrt(temperature/(2*damp*dt)); 
    rot_bd_prefactor= sqrt(temperature/(2*damp*dt)); 
    
    /********for FENE springs*********/ 
    max_ext     = max_ext_ratio*mlen;
    eps_ext     = 0.01*max_ext;
    /********************************/

    shear       = 0;
    tension     = 0;
    force       = {{0,0}}; // force on the spring  
    b_eng       = {{0,0}}; // filament / xlink bending energy
    
    b_force[0]  = {{0,0}}; //b_force[0] = bending force on head 0 due to h0-h1-spring angle in cartesian coords
    b_force[1]  = {{0,0}}; //b_force[1] = bending force on head 1 due to h1-h0-spring angle in cartesian coords
    
    kinetic_energy = 0;
    pos_rat = {{0, 0}}; // pos_rat = distance from pointy end / segment length-- by default 0
                        // i.e., if l_index[hd] = j, then pos_rat[hd]*|disp[hd]| is the distance to the "j+1"th bead
    
    array<double, 2> posH0 = boundary_check(0, pos[0], pos[1]); 
    array<double, 2> posH1 = boundary_check(1, pos[0]+pos[2], pos[1]+pos[3]); 
    hx[0] = posH0[0];
    hy[0] = posH0[1];
    hx[1] = posH1[0];
    hy[1] = posH1[1];
    
    if (state[0] == 1){
        update_pos_rat(0);
        ldir_bind[0] = filament_network->get_direction(f_index[0], l_index[0]);

    }
    if (state[1] == 1){
        update_pos_rat(1);
        ldir_bind[1] = filament_network->get_direction(f_index[1], l_index[1]);
    }
    
    this->update_angle();
    this->update_th_from_pos();
    this->update_direc();
    cm = pos_bc(BC, filament_network->get_delrx(), dt, fov, {{0, 0}}, {{hx[0]+0.5*mld*direc[0], hy[0]+0.5*mld*direc[1]}});
   
    //force can be non-zero and angle is determined from disp vector
    this->update_force();
    
    ldir_bind[0] = {{0,0}};
    ldir_bind[1] = {{0,0}};
    bind_disp[0] = {{0,0}};
    bind_disp[1] = {{0,0}};
    
    bind_rot = {{0,0}};

    at_barbed_end = {{false, false}};


    prv_rnd_x = 0;
    prv_rnd_y = 0;
    prv_rnd_th = 0;

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

        //project bending force onto filaments
        this->filament_update_hd(0, {{ b_force[0][0], b_force[0][1] }});
        this->filament_update_hd(1, {{ b_force[1][0], b_force[1][1] }});
        
        //reset bending force
        b_force[0] = {{0,0}};
        b_force[1] = {{0,0}};
        //force      = {{0,0}};
       
        //enforce constraints
        update_shake_force();
    }
    else
    {
        b_force[0] = {{0,0}};
        b_force[1] = {{0,0}};
    }

    //update_angle(); // need to recalculate this, since heads might've moved (this doesn't make sense)

}

//Measure distance to FARTHER END of bead filament that the spacer is bound to
int spacer::get_further_end(int hd, int findex, int lindex)
{
    return (pos_rat[hd] > 0.5);
}

array<double, 2> spacer::disp_from_bead(int hd, int findex, int aindex)
{
  return rij_bc(BC,
          filament_network->get_filament(findex)->get_bead(aindex)->get_xcm() - hx[hd],
          filament_network->get_filament(findex)->get_bead(aindex)->get_ycm() - hy[hd],
          fov[0], fov[1], filament_network->get_delrx());
}

void spacer::update_bending(int hd)
{
    array<double, 2> delr1, delr2;
    double f1[2], f3[2];
    double rsq1,rsq2,r1,r2,c,s,dth,a,a11,a12,a22;

    int bead_further_end = get_further_end(hd, f_index[hd], l_index[hd]);

    // 1st bond
    delr1 = disp_from_bead(hd, f_index[hd], l_index[hd] + bead_further_end); 
    rsq1  = delr1[0]*delr1[0] + delr1[1]*delr1[1];
    r1    = sqrt(rsq1);

    // 2nd bond
    delr2 = {{pow(-1, hd)*disp[0], pow(-1, hd)*disp[1]}};
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

    filament_network->update_forces(f_index[hd], l_index[hd] + bead_further_end, f1[0], f1[1]);
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

void spacer::brownian_relax(int hd){
    
    if ( hd == 0 ){
        if ( state[1] == 0 )
            this->update_unattached();
        else 
            this->update_one_attached(1);
    }
    else
    {
        if (state[0] == 1)
            this->update_one_attached();
        //else do nothing because already updated unnattached
    }
    
    update_hx_hy();

}
void spacer::update_unattached()
{
    //cout<<"\nDEBUG: using spacer brownian_relax()";
    
    double new_rnd_x = rng_n(), new_rnd_y = rng_n(), new_rnd_th = rng_n();
    double vx =  bd_prefactor*(new_rnd_x + prv_rnd_x);
    double vy =  bd_prefactor*(new_rnd_y + prv_rnd_y);
    double w  =  rot_bd_prefactor*(new_rnd_th + prv_rnd_th);
    
    //kinetic_energy = vx*vx + vy*vy;    
    cm = pos_bc(BC, filament_network->get_delrx(), dt, fov, {{vx, vy}}, {{cm[0]+vx*dt, cm[1]+vy*dt}});
    th = th + w*dt;
    update_direc();

    prv_rnd_x = new_rnd_x;
    prv_rnd_y = new_rnd_y;
    prv_rnd_th = new_rnd_th;

}

void spacer::update_one_attached(int at_hd)
{
    double new_rnd_th = rng_n();
    double w =  rot_bd_prefactor*(new_rnd_th + prv_rnd_th);
    
    th = th + w*dt/2;
    update_direc();
   
    int hd_sgn = pow(-1, at_hd);
    array<double, 2> dr = {{0.5*hd_sgn*mld*direc[0], 0.5*hd_sgn*mld*direc[1]}};
    array<double, 2> cm_new = {{hx[at_hd] + dr[0], hy[at_hd] + dr[1]}};
    cm = pos_bc(BC, filament_network->get_delrx(), dt, fov, {{0, 0}}, cm_new);
    
//    array<double, 2> un_hd_pos = {{cm[0] + dr[0], cm[1] + dr[1]}};
//    un_hd_pos = pos_bc(BC, filament_network->get_delrx(), dt, fov, {{0,0}}, un_hd_pos);

    prv_rnd_th = new_rnd_th;
}

void spacer::update_direc()
{
    direc = {{cos(th), sin(th)}}; 
}

void spacer::update_th_from_pos()
{
    array<double, 2> disp = rij_bc(BC, hx[1]-hx[0], hy[1]-hy[0], fov[0], fov[1], filament_network->get_delrx()); 
    th = atan2(disp[1], disp[0]);
}

void spacer::update_hx_hy()
{
    array<double, 2> dr = {{0.5*mld*direc[0], 0.5*mld*direc[1]}}; 
    array<double, 2> r0 = pos_bc(BC, filament_network->get_delrx(), dt, fov, {{0,0}}, {{cm[0] - dr[0], cm[1]-dr[1]}});
    array<double, 2> r1 = pos_bc(BC, filament_network->get_delrx(), dt, fov, {{0,0}}, {{cm[0] + dr[0], cm[1]+dr[1]}});
    hx = {{r0[0], r1[0]}};
    hy = {{r0[1], r1[1]}};
}

//check for attachment of unbound heads given head index (0 for head 1, and 1 for head 2)
bool spacer::attach(int hd)
{
    double not_off_prob = 0, onrate = kon;
    double mf_rand = rng(0,1.0);
    array<double, 2> intPoint;
    if (state[pr(hd)] == 1)
        onrate = kon2;
   
//    cout<<"\nDEBUG: head = "<<hd;
    set<pair<array<double, 2>, array<int, 2> > > dist_sorted = filament_network->get_binding_points({{hx[pr(hd)], hy[pr(hd)]}}, mld);

    if(!dist_sorted.empty()){
        

        for (set<pair<array<double,2>, array<int, 2>>>::iterator it=dist_sorted.begin(); it!=dist_sorted.end(); ++it)
        {
            //head can't bind to the same filament link the other head is bound to
            //cout<<"\nDEBUG: trying to bind head "<<hd<<" to filament "<<it->second.at(0);
            if(allowed_bind(hd, it->second)){
                
                intPoint = it->first;
                not_off_prob += metropolis_prob(hd, it->second, intPoint, onrate);
                //cout<<"\nDEBUG: probability of binding head "<<hd<<" to filament "<<it->second.at(0) << " = "<<not_off_prob;
                         
                if (mf_rand < not_off_prob) 
                {
                    //cout<<"\nDEBUG: bound head "<<hd<<" to filament "<<it->second.at(0);
                    //update state
                    state[hd] = 1;
                    f_index[hd] = (it->second).at(0);
                    this->set_l_index(hd, (it->second).at(1));
                    
                    //update head position
                    hx[hd] = intPoint[0];
                    hy[hd] = intPoint[1];
                    double thprev = th;
                    update_th_from_pos();
                    update_direc();

                    //record rotation of head for future unbinding move
                    bind_rot[hd] = th - thprev;

                    //update relative head position
                    update_pos_rat(hd);
                    
                    //(even if its at the barbed end upon binding, could have negative velocity, so always set this to false, until it steps)
                    at_barbed_end[hd] = false; 
                    
                    if (state[pr(hd)] == 1)
                        disp_prev = disp;

                    return true;
                }
            }
        }
    }	
    return false;
} 


void spacer::filament_update()
{
    if (state[0] == 1 && state[1] == 1){

        this->filament_update_hd(0, {{ force[0]/* + b_force[0][0]*/,  force[1] /*+ b_force[0][1]*/}});
        this->filament_update_hd(1, {{-force[0]/* + b_force[1][0]*/, -force[1] /*+ b_force[1][1]*/}});

        //reset bending force
        //b_force[0] = {{0,0}};
        //b_force[1] = {{0,0}};
        force      = {{0,0}};
    }
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
//    double stretch  = dist_bc(BC, newpos[0] - hx[pr(hd)], newpos[1] - hy[pr(hd)], fov[0], fov[1], filament_network->get_delrx()) - mld; 
    
    array<double, 2> delr1, delr2;
    double r1, r2, c, dth, bend_eng = 0;
    
    if (state[hd] == 0 && state[pr(hd)] == 1) { //it's trying to attach
        
//        delr1 = disp_from_bead(hd, it->second.at(0), it->second.at(1) + get_further_end(hd, it->second.at(0), it->second.at(1))); 
        delr1 = disp_from_bead(hd, fl_idx[0], fl_idx[1] + get_further_end(hd, fl_idx[0], fl_idx[1])); 
        r1  = sqrt(delr1[0]*delr1[0] + delr1[1]*delr1[1]);

        // 2nd bond
        delr2 = {{pow(-1, hd)*disp[0], pow(-1, hd)*disp[1]}};
        r2  = sqrt(delr2[0]*delr2[0] + delr2[1]*delr2[1]);

        // cos
        c = (delr1[0]*delr2[0] + delr1[1]*delr2[1]) / (r1*r2);

        if (c > 1.0) c = 1.0;
        if (c < -1.0) c = -1.0;

        dth = acos(c) - th0;
        bend_eng = kb*dth*dth/(r1+r2);
        //cout<<"\nDEBUG: bend_eng = "<<bend_eng; 
    }
    
//    double delE = 0.5*mk*stretch*stretch + bend_eng - this->get_stretching_energy() - b_eng[hd];
    double delE = bend_eng - b_eng[hd];
    if( delE > 0 )
        prob *= exp(-delE/temperature);
    
    return prob;
}

bool spacer::allowed_bind(int hd, array<int, 2> fl_idx)
{
//    cout<<"\nDEBUG: using spacer allowed bind";
    return (fl_idx[0] != f_index[pr(hd)]);
}

array<double, 2> spacer::generate_off_pos(int hd)
{
    
    double thnew = th + bind_rot[hd];
    return pos_bc(BC, filament_network->get_delrx(), 0, fov, {{0,0}}, {{hx[pr(hd)] + mld*cos(thnew), hy[pr(hd)] + mld*sin(thnew)}}); 

} 

/* ---------------------------------------------------------------------- */

void spacer::update_shake_force()
{
    // r01 = distance vec between atoms, with PBC
    array<double, 2> r01 = disp;

    //get unconstrained positions of beads of links that xlink is connected to
    //(the convention for pos_rat is measured from the pointed to the barbed end, hence pos00 is closer to pointed end)
    array<double, 2> pos00 = filament_network->get_filament(f_index[0])->get_predicted_position(l_index[0]+1);
    array<double, 2> pos01 = filament_network->get_filament(f_index[0])->get_predicted_position(l_index[0]);
    array<double, 2> pos10 = filament_network->get_filament(f_index[1])->get_predicted_position(l_index[1]+1);
    array<double, 2> pos11 = filament_network->get_filament(f_index[1])->get_predicted_position(l_index[1]);

    array<double, 2> disp0 = rij_bc(BC, pos01[0]-pos00[0], pos01[1]-pos00[1], fov[0], fov[1], filament_network->get_delrx());
    array<double, 2> disp1 = rij_bc(BC, pos11[0]-pos10[0], pos11[1]-pos10[1], fov[0], fov[1], filament_network->get_delrx());
    
    //calculate unconstrained position of crosslinker head using unconstrained positions of filament beads 
    array<double, 2> unconstrained_h0 = pos_bc(BC, filament_network->get_delrx(), dt, fov, {{0, 0}}, 
            {{pos00[0] + pos_rat[0]*disp0[0], pos00[1] + pos_rat[0]*disp0[1]}});
    array<double, 2> unconstrained_h1 = pos_bc(BC, filament_network->get_delrx(), dt, fov, {{0, 0}}, 
            {{pos10[0] + pos_rat[1]*disp1[0], pos10[1] + pos_rat[1]*disp1[1]}});
    
    //enforce constraint on crosslinker 
    // s01 = distance vec after unconstrained update, with PBC
    array<double, 2> s01 = rij_bc(BC, 
            unconstrained_h1[0] - unconstrained_h0[0], 
            unconstrained_h1[1] - unconstrained_h0[1], fov[0], fov[1], filament_network->get_delrx()); 


    // scalar distances between atoms
    double r01sq = r01[0]*r01[0] + r01[1]*r01[1]; //+ r01[2]*r01[2];
    double s01sq = s01[0]*s01[0] + s01[1]*s01[1]; //+ s01[2]*s01[2];

    double invm = 1/(damp*dt);
    double a = invm*invm * r01sq; //not sure why this is r01sq instead of mld^2 like in Frenk&Smit p.413
    double b = 2.0 * invm * dot(s01, r01);
    double c = s01sq - mld*mld;

    // error check
    double determ = b*b - 4.0*a*c;
    if (determ < 0.0) {
        cout<<"\nWARNING: Shake determinant < 0.0; setting to 0";
        determ = 0.0;
    }

    // exact quadratic solution for lamda

    double lamda,lamda1,lamda2;
    lamda1 = (-b+sqrt(determ)) / (2.0*a);
    lamda2 = (-b-sqrt(determ)) / (2.0*a);

    if (fabs(lamda1) <= fabs(lamda2)) 
        lamda = lamda1;
    else 
        lamda = lamda2;

    // update forces if atom is owned by this processor
    lamda /= (dt*dt);

    //update_force
    force[0] = -lamda*r01[0];
    force[1] = -lamda*r01[1];

    //project force onto filaments
}
array<double, 2> spacer::get_bending_energy()
{
    return b_eng;
}
