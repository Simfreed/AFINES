/*------------------------------------------------------------------
 filament.cpp : object describing a worm-like chain filament
 
 Copyright (C) 2016 
 Created by: Simon Freedman, Shiladitya Banerjee, Glen Hocky, Aaron Dinner
 Contact: dinner@uchicago.edu

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version. See ../LICENSE for details. 
-------------------------------------------------------------------*/

#include "filament.h"
#include "bead.h"
#include "globals.h"
//using namespace std;
filament::filament(){

    fov[0] = 50;
    fov[1] = 50;
    nq[0] = 100;
    nq[1] = 100;
    dt = 0.001;
    temperature = 0;
    fracture_force = 1000000;
    BC = "REFLECTIVE";
    kinetic_energy = 0;
    gamma = 0;
    delrx=0;
    damp = infty;
    y_thresh=2;
    bd_prefactor = sqrt(temperature/(2*dt*damp));
}

filament::filament(array<double, 2> myfov, array<int, 2> mynq, double deltat, double temp, double shear, 
        double frac, double bending_stiffness, string bndcnd)
{
    fov             = myfov;
    nq              = mynq;
    dt              = deltat;
    temperature     = temp;
    gamma           = shear;
    delrx           = 0;
    fracture_force  = frac;
    kb              = bending_stiffness;
    BC              = bndcnd;
    kinetic_energy  = 0;
    damp            = infty;
    y_thresh        = 1;
    bd_prefactor = sqrt(temperature/(2*dt*damp));

}

filament::filament(array<double, 3> startpos, int nbead, array<double, 2> myfov, array<int, 2> mynq, double visc, 
        double deltat, double temp, bool isStraight, double beadRadius, double spring_length, double stretching_stiffness,
        double max_ext_ratio, double bending_stiffness, double frac_force, string bdcnd)
{
    
    fov = myfov;
    nq = mynq;
    dt = deltat;
    temperature = temp;
    gamma = 0;
    delrx = 0;
    fracture_force = frac_force;
    BC = bdcnd;
    kb = bending_stiffness;
    kinetic_energy  = 0;
    
    damp = 6*pi*beadRadius*visc;
    y_thresh = 1;
    
    bd_prefactor = sqrt(temperature/(2*dt*damp));

    double phi, variance;
    array<double, 2> next_pos;
    //the start of the polymer: 
    beads.push_back(new bead( startpos[0], startpos[1], beadRadius, visc));
    prv_rnds.push_back({0,0});
    phi = startpos[2];
    
    if (temp != 0) variance = temp/bending_stiffness;
    else variance = 0;

    for (int j = 1; j < nbead; j++) {

        //next_pos = boundary_check(j-1, beads[j-1]->get_xcm() + spring_length*cos(phi), beads[j-1]->get_ycm() + spring_length*sin(phi));
        next_pos = pos_bc(BC, delrx, dt, fov, 
                {spring_length*cos(phi)/dt, spring_length*sin(phi)/dt},
                {beads[j-1]->get_xcm() + spring_length*cos(phi), beads[j-1]->get_ycm() + spring_length*sin(phi)});
        beads.push_back( new bead(next_pos[0], next_pos[1], beadRadius, visc) );
        prv_rnds.push_back({0,0});
        springs.push_back( new spring(spring_length, stretching_stiffness, max_ext_ratio, this, {j-1, j}, fov, nq) );  
        springs[j-1]->step(BC, delrx);  
        springs[j-1]->update_force(BC, delrx);
        
        // Calculate the Next angle on the bead polymer
        if (!isStraight) phi += rng_n(0, variance);
    
    }
   
}

filament::filament(vector<bead *> beadvec, array<double, 2> myfov, array<int, 2> mynq, double spring_length, 
        double stretching_stiffness, double max_ext_ratio, double bending_stiffness, 
        double deltat, double temp, double frac_force, double g, string bdcnd)
{
    
    dt = deltat;
    temperature = temp;
    fracture_force = frac_force;
    gamma = g; 
    delrx = 0;
    BC = bdcnd;
    kb = bending_stiffness;
    fov = myfov;
    nq = mynq;
    y_thresh = 1;
    kinetic_energy = 0;

    if (beadvec.size() > 0)
    {
        beads.push_back(new bead(*(beadvec[0])));
        prv_rnds.push_back({0,0});
        damp = beads[0]->get_friction();
    }
    
    //spring em up
    if (beadvec.size() > 1){
        for (unsigned int j = 1; j < beadvec.size(); j++) {

            beads.push_back(new bead(*(beadvec[j])));
            springs.push_back( new spring(spring_length, stretching_stiffness, max_ext_ratio, this, {(int)j-1, (int)j}, fov, nq) );  
            springs[j-1]->step(BC, delrx);
            springs[j-1]->update_force(BC, delrx);
            prv_rnds.push_back({0,0});
            
        }
    }
    
    bd_prefactor = sqrt(temperature/(2*dt*damp));

}

filament::~filament(){
    
    //cout<<"DELETING FILAMENT\n";
    int nr = beads.size(), nl = springs.size();
    for (int i = 0; i < nr; i ++)
    {    
        //cout<<"\nDEBUG: deleting pointer "<<beads[i];     
        delete beads[i];
    }
    for (int i = 0; i < nl; i ++)
        delete springs[i];
    
    beads.clear();
    springs.clear();
    prv_rnds.clear();
}

void filament::add_bead(bead * a, double spring_length, double stretching_stiffness, double max_ext_ratio){
    
    beads.push_back(new bead(*a));
    prv_rnds.push_back({0,0});    
    if (beads.size() > 1){
        int j = (int) beads.size() - 1;
        springs.push_back( new spring(spring_length, stretching_stiffness, max_ext_ratio, this, {j-1,  j}, fov, nq ) );  
        springs[j-1]->step(BC, delrx);
    }
    if (damp == infty)
        damp = a->get_friction();
}

vector<vector<array<int,2> > > filament::get_quadrants()
{
    //should return a map between bead and x, y coords of quadrant
    vector<vector<array<int,2> > > quads;
    for (unsigned int i=0; i < springs.size(); i++){ 
        springs[i]->quad_update(BC, delrx);
        quads.push_back(springs[i]->get_quadrants());
    }
    
    return quads;
}
/*
multimap<int, array<int, 2> > filament::get_quadrants()
{
    //should return a map between bead and x, y coords of quadrant
    multipmap<int, array<int, 2> > quads;
    for (unsigned int i=0; i < springs.size(); i++){ 
        springs[i]->quad_update(BC, delrx);
        quads.push_back(springs[i]->get_quadrants());
    }
    
    return quads;
}*/

void filament::set_y_thresh(double y){
    y_thresh = y;
}

void filament::update_positions()
{
    double vx, vy;
    array<double, 2> new_rnds;
    array<double, 2> newpos;
    kinetic_energy = 0;  
    double top_y = y_thresh*fov[1]/2.; 
    int sa = int(beads.size());
    int la = int(springs.size());
    for (int i = 0; i < sa; i++){

        if (fabs(beads[i]->get_ycm()) > top_y) continue;
     
        new_rnds = {rng_n(), rng_n()};
        vx  = (beads[i]->get_force()[0])/damp  + bd_prefactor*(new_rnds[0] + prv_rnds[i][0]);
        vy  = (beads[i]->get_force()[1])/damp  + bd_prefactor*(new_rnds[1] + prv_rnds[i][1]);
//        cout<<"\nDEBUG: Fx("<<i<<") = "<<beads[i]->get_force()[0]<<"; v = ("<<vx<<" , "<<vy<<")";
       
        prv_rnds[i] = new_rnds;
        //cout<<"\nDEBUG: bead force = ("<<beads[i]->get_force()[0]<<" , "<<beads[i]->get_force()[1]<<")";
        kinetic_energy += vx*vx + vy*vy;
        newpos = pos_bc(BC, delrx, dt, fov, {vx, vy}, {beads[i]->get_xcm() + vx*dt, beads[i]->get_ycm() + vy*dt});
        beads[i]->set_xcm(newpos[0]);
        beads[i]->set_ycm(newpos[1]);
        beads[i]->reset_force(); 
    }

    for (int i = 0; i < la; i++)
        springs[i]->step(BC, delrx);

}

void filament::update_positions_range(int lo, int hi)
{
    double vx, vy;
    array<double, 2> new_rnds;
    array<double, 2> newpos;
    kinetic_energy = 0;  
    double top_y = y_thresh*fov[1]/2.; 

    int low = max(0, lo);
    int high = min(hi, (int)beads.size());

    for (int i = low; i < high; i++){
       
        if (fabs(beads[i]->get_ycm()) > top_y) continue;
     
        new_rnds = {rng_n(), rng_n()};
        vx  = (beads[i]->get_force()[0])/damp  + bd_prefactor*(new_rnds[0] + prv_rnds[i][0]);
        vy  = (beads[i]->get_force()[1])/damp  + bd_prefactor*(new_rnds[1] + prv_rnds[i][1]);
//        cout<<"\nDEBUG: Fx("<<i<<") = "<<beads[i]->get_force()[0]<<"; v = ("<<vx<<" , "<<vy<<")";
       
        prv_rnds[i] = new_rnds;
        //cout<<"\nDEBUG: bead force = ("<<beads[i]->get_force()[0]<<" , "<<beads[i]->get_force()[1]<<")";
        kinetic_energy += vx*vx + vy*vy;
        //newpos = boundary_check(i, beads[i]->get_xcm() + vx*dt, beads[i]->get_ycm() + vy*dt); 
        newpos = pos_bc(BC, delrx, dt, fov, {vx, vy}, {beads[i]->get_xcm() + vx*dt, beads[i]->get_ycm() + vy*dt});
        beads[i]->set_xcm(newpos[0]);
        beads[i]->set_ycm(newpos[1]);
        beads[i]->reset_force(); 
    }

    for (unsigned int i = 0; i < springs.size(); i++)
        springs[i]->step(BC, delrx);

}

vector<filament *> filament::update_stretching(double t)
{
    vector<filament *> newfilaments;

    if(springs.size() == 0)
        return newfilaments;
   
    for (unsigned int i=0; i < springs.size(); i++) {
        springs[i]->update_force(BC, delrx);
        //springs[i]->update_force_fraenkel_fene(BC, delrx);
        if (hypot(springs[i]->get_force()[0], springs[i]->get_force()[1]) > fracture_force){
            newfilaments = this->fracture(i);
            break;
        }
        else 
            springs[i]->filament_update();
    }
    
    return newfilaments;
}

bead * filament::get_bead(int i)
{
    try
    {
        return beads[i];
    }
    catch (int e)
    {
        cout<<"\nDEBUG: an exception occured while returning the beads[ "<<i<<"]";
        bead * a;
        return a;
    }
}

spring * filament::get_spring(int i)
{
    return springs[i];
}

void filament::update_delrx(double shear_dist){
    delrx = shear_dist;
}

void filament::update_shear(double t){
    
    double local_shear;
    for (unsigned int i = 0; i < beads.size(); i++){
        local_shear = delrx * beads[i]->get_ycm() / fov[1];
        beads[i]->set_xcm(beads[i]->get_xcm() + local_shear);
        //cout<<"\nDEBUG: local_shear = "<<local_shear;
    }
}

void filament::update_d_strain(double g){
    
    for (unsigned int i = 0; i < beads.size(); i++){
        beads[i]->set_xcm(beads[i]->get_xcm() + g * beads[i]->get_ycm() / fov[1]);
    }
}

void filament::update_forces(int index, double f1, double f2)
{
    beads[index]->update_force(f1,f2);
}

void filament::pull_on_ends(double f)
{
    if (beads.size() < 2) return;
    int last = beads.size() - 1; 
    array<double, 2> dr = rij_bc(BC, beads[last]->get_xcm() - beads[0]->get_xcm(),  
                                     beads[last]->get_ycm() - beads[0]->get_ycm(), fov[0], fov[1], delrx);
    double ang = atan2( dr[1], dr[0]);
    
    beads[ 0  ]->update_force(-0.5*f*cos(ang), -0.5*f*sin(ang));
    beads[last]->update_force( 0.5*f*cos(ang),  0.5*f*sin(ang));
}

void filament::affine_pull(double f)
{
    if (beads.size() < 2) return;
    int last = beads.size() - 1; 
    array<double, 2> dr = rij_bc(BC, beads[last]->get_xcm() - beads[0]->get_xcm(),  
                                     beads[last]->get_ycm() - beads[0]->get_ycm(), fov[0], fov[1], delrx);
    double ang = atan2( dr[1], dr[0]);
    //cout<<"\nDEBUG: angle = "<<ang;
    double frac, fcos = f*cos(ang), fsin = f*sin(ang);

    for (int i = 0; i <= last; i++){
        frac = (double(i)/double(last)-0.5);
        beads[i]->update_force(frac*fcos, frac*fsin);
    }
}

void filament::set_shear(double g){
    gamma = g;
    max_shear = gamma*fov[1]*0.5;
}

string filament::write_beads(int fil){
    string all_beads;
    for (unsigned int i =0; i < beads.size(); i++)
    {
        all_beads += beads[i]->write() + "\t" + std::to_string(fil);
    }

    return all_beads;
}

string filament::write_springs(int fil){
    string all_springs;
    for (unsigned int i =0; i < springs.size(); i++)
    {
        all_springs += springs[i]->write(BC, delrx) + "\t" + std::to_string(fil);
    }

    return all_springs;
}

string filament::write_thermo(int fil)
{
    return "\n" + std::to_string(this->get_kinetic_energy()) + \
        "\t" + std::to_string(this->get_potential_energy()) + \
        "\t" + std::to_string(this->get_total_energy()) + "\t" + std::to_string(fil);
}

vector<bead *> filament::get_beads(unsigned int first, unsigned int last)
{
    vector<bead *> newbeads;
    for (unsigned int i = first; i < last; i++)
    {
        if (i >= beads.size())
            break;
        else
            newbeads.push_back(new bead(*(beads[i])));
    }
    return newbeads;
}

vector<filament *> filament::fracture(int node){

    vector<filament *> newfilaments;
    cout<<"\n\tDEBUG: fracturing at node "<<node;
    
    if(springs.size() == 0)
        return newfilaments;

    vector<bead *> lower_half = this->get_beads(0, node+1);
    vector<bead *> upper_half = this->get_beads(node+1, beads.size());

    if (lower_half.size() > 0)
        newfilaments.push_back(
                new filament(lower_half, fov, nq, springs[0]->get_l0(), springs[0]->get_kl(), springs[0]->get_fene_ext(), kb, 
                    dt, temperature, fracture_force, gamma, BC));
    if (upper_half.size() > 0)
        newfilaments.push_back(
                new filament(upper_half, fov, nq, springs[0]->get_l0(), springs[0]->get_kl(), springs[0]->get_fene_ext(), kb, 
                    dt, temperature, fracture_force, gamma, BC));

    for (int i = 0; i < (int)(lower_half.size()); i++) delete lower_half[i];
    for (int i = 0; i < (int)(upper_half.size()); i++) delete upper_half[i];
    
    lower_half.clear();
    upper_half.clear();
    
    return newfilaments;

}

bool filament::operator==(const filament& that){
    
    if (beads.size() != that.beads.size() || springs.size() != that.springs.size())
        return false;

    for (unsigned int i = 0; i < beads.size(); i++)
        if (!(*(beads[i]) == *(that.beads[i])))
            return false;
    
    for (unsigned int i = 0; i < springs.size(); i++)
        if (!(springs[i]->is_similar(*(that.springs[i]))))
            return false;

    return (this->fov[0] == that.fov[0] && this->fov[1] == that.fov[1] && 
            this->nq[0] == that.nq[0] && this->nq[1] == that.nq[1] &&
            this->gamma == that.gamma && this->temperature == that.temperature &&
            this->dt == that.dt && this->fracture_force == that.fracture_force);

}

string filament::to_string(){
    
    // Note: not including springs in to_string, because spring's to_string includes filament's to_string
    char buffer[200];
    string out = "";

    for (unsigned int i = 0; i < beads.size(); i++)
        out += beads[i]->to_string();

    sprintf(buffer, "fov = (%f, %f)\tnq = (%d, %d)\tgamma = %f\ttemperature = %f\tdt = %f\tfracture_force=%f\n",
            fov[0], fov[1], nq[0], nq[1], gamma, temperature, dt, fracture_force);
   
    return out + buffer; 

}

string filament::get_BC(){
    return BC; 
}

void filament::set_BC(string s){
    this->BC = s;
}

inline double filament::angle_between_springs(int i, int j){
    double theta = springs[i]->get_angle() - springs[j]->get_angle();
    return theta - 2*pi*floor(theta/(2*pi)+0.5); // Keep angles between -Pi and Pi
}

/* --------------------------------------------------------------------------------- */
/* copied and pasted and edited the following code from LAMMPS src/angle_harmonic.cpp*/
/* --------------------------------------------------------------------------------- */
void filament::lammps_bending_update()
{
  array<double, 2> delr1, delr2;
  double f1[2], f3[2];
  double rsq1,rsq2,r1,r2,c,s,a,a11,a12,a22;
    
  // 1st bond
  delr1 = springs[0]->get_neg_disp();
  rsq1  = delr1[0]*delr1[0] + delr1[1]*delr1[1];
  r1    = sqrt(rsq1);

  for (int n = 0; n < int(springs.size())-1; n++) {
    
    // 2nd bond
    delr2 = springs[n+1]->get_disp();
    rsq2  = delr2[0]*delr2[0] + delr2[1]*delr2[1];
    r2    = sqrt(rsq2);

    // angle (cos and sin)
    c = (delr1[0]*delr2[0] + delr1[1]*delr2[1]) / (r1*r2);
    
    if (c > 1.0) c = 1.0;
    if (c < -1.0) c = -1.0;

    s = sqrt(1.0 - c*c);
    if (s < maxSmallAngle) s = maxSmallAngle;

    // force
    a   = -kb * (acos(c) - pi) / s; //Note, in this implementation, Lp = kb/kT 
    a11 = a*c / rsq1;
    a12 = -a / (r1*r2);
    a22 = a*c / rsq2;
  
    f1[0] = a11*delr1[0] + a12*delr2[0];
    f1[1] = a11*delr1[1] + a12*delr2[1];
    f3[0] = a22*delr2[0] + a12*delr1[0];
    f3[1] = a22*delr2[1] + a12*delr1[1];
    //cout<<"\nDEBUG: f1x, f1y, f3x f3y = "<<f1[0]<<" , "<<f1[1]<<" , "<<f3[0]<<" , "<<f3[1];
                
    // apply force to each of 3 atoms
    beads[n  ]->update_force(f1[0], f1[1]);
    beads[n+1]->update_force(-f1[0] - f3[0], -f1[1] - f3[1]);
    beads[n+2]->update_force(f3[0], f3[1]);

    // 1st bond, next iteration
    delr1 = {-delr2[0], -delr2[1]};
    rsq1  = rsq2;
    r1    = r2;

  }
}

void filament::fwd_bwd_bending_update()
{
  array<double, 2> delr1, delr2;
  double rsq1,rsq2,r1,r2,c,c1,c2,c3,cost,sint;//,s,coeff;

  // 1st bond
  delr1 = springs[0]->get_disp();
  rsq1  = delr1[0]*delr1[0] + delr1[1]*delr1[1];///Caa
  r1    = sqrt(rsq1);

  for (int n = 1; n < int(beads.size())-1; n++) {

    // 2nd bond
    delr2 = springs[n]->get_disp();
    rsq2  = delr2[0]*delr2[0] + delr2[1]*delr2[1];//Ca+1a+1
    r2    = sqrt(rsq2); 

    // angle (cos and sin)

    if (r1 < eps || r2 < eps ) // bad news, probably inaccurate results
        continue;

    c  = delr1[0]*delr2[0] + delr1[1]*delr2[1];//Caa+1
    
    cost = c /(r1*r2); // negative because cos(pi-theta) = -cos(theta)

    if (cost > 1.0) cost = 1.0;
    if (cost < -1.0) cost = -1.0;

    sint = sqrt(1.0 - cost*cost);
    if (sint < maxSmallAngle) sint = maxSmallAngle;

    c1 = kb*acos(cost)/(sint*r1*r2);
    c2 = c/rsq2;
    c3 = c/rsq1;
    beads[n]->update_force(
            c1*( delr2[0]*(1+c2) - delr1[0]*(1+c3) ),
            c1*( delr2[1]*(1+c2) - delr1[1]*(1+c3) )
            );

    /*
    
    // force
    coeff = -2*kb*acos(cost)/s;

    // Allen and Tildesley Appendix C, eq C.13 middle
    // (top and bottom of C.13 cancels out when going fwd and bwd; middle doubles

    beads[n]->update_force(
            coeff*( c*( delr2[0] / rsq2 - delr1[0] / rsq1 ) - ( delr2[0] - delr1[0] )/( r1*r2 ) ),
            coeff*( c*( delr2[1] / rsq2 - delr1[1] / rsq1 ) - ( delr2[1] - delr1[1] )/( r1*r2 ) )
            );
    */    
    // 1st bond of next iteration
    delr1 = delr2;
    rsq1  = rsq2;
    r1    = r2;
  }
}

//wrapper, for fwd_bending_update (and bwd bending update if I ever make it)
void filament::update_bending(double t)
{
    if(springs.size() > 1 && kb > 0){
          this->lammps_bending_update();
//          this->fwd_bwd_bending_update();
    }
}

int filament::get_nbeads(){
    return beads.size();
}

int filament::get_nsprings(){
    return springs.size();
}

double filament::get_bending_energy(){
    
    double sum = 0, theta;

    if (springs.size() < 2) return 0;

    for (unsigned int i = 0; i < springs.size() - 1; i++)
    {
        theta = angle_between_springs(i+1, i);//springs[i+1]->get_angle() - springs[i]->get_angle();
        sum += theta*theta;
    }
    
    return kb*sum/2;

}

double filament::get_stretching_energy()
{
    
    double u = 0;
    for (unsigned int i = 0; i < springs.size(); i++)
        //u += springs[i]->get_stretching_energy_fene(BC, delrx);
        u += springs[i]->get_stretching_energy();
    
    return u;

}

double filament::get_kinetic_energy()
{
    return kinetic_energy;
}

double filament::get_potential_energy()
{
    return this->get_stretching_energy() + this->get_bending_energy();
}

double filament::get_total_energy()
{
    return this->get_potential_energy() + this->get_kinetic_energy();
}

array<double, 2> filament::get_bead_position(int n)
{
    return {beads[n]->get_xcm(), beads[n]->get_ycm()};
}

void filament::print_thermo()
{
    cout<<"\tKE = "<<this->get_kinetic_energy()<<"\tPE = "<<this->get_potential_energy()<<\
        "\tTE = "<<this->get_total_energy();
}

double filament::get_end2end()
{
    if (beads.size() < 2) 
        return 0;
    else 
        return dist_bc(BC, beads[beads.size() - 1]->get_xcm() - beads[0]->get_xcm(),  
                           beads[beads.size() - 1]->get_ycm() - beads[0]->get_ycm(), fov[0], fov[1], delrx);
} 
