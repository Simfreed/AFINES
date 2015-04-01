/*
 *  actin.cpp
 *  
 *
 *  Created by Shiladitya Banerjee on 9/3/13.
 *  Modified by Simon Freedman 9/2014
 *  Copyright 2013 University of Chicago. All rights reserved.
 *
 */

#include "globals.h"
#include "link_ensemble.h"
#include "Link.h"
#include "filament_ensemble.h"
//actin network class

template <class filament_type> 
filament_ensemble<filament_type>::filament_ensemble(){}

template <class filament_type> 
filament_ensemble<filament_type>::~filament_ensemble(){ 
    cout<<"DELETING FILAMENT_ENSEMBLE\n";
    
    int s = network.size();
    for (int i = 0; i < s; i++){
        delete network[i];
    }
    
    //network.clear();
};

template <class filament_type> 
void filament_ensemble<filament_type>::quad_update()
{
    quad_fils.clear();
    vector< vector< array<int, 2> > > filament_quads;

    for (unsigned int i=0; i < network.size(); i++) { //Loop over filaments
        
        filament_quads = network[i]->get_quadrants(); 
        
        for (unsigned int j=0; j < filament_quads.size(); j++){ //Loop over links
            
            for (unsigned int k = 0; k  < filament_quads[j].size(); k++){ //Loop over quadrants of a Link   
                
                quad_fils[ filament_quads[j][k] ].push_back({(int)i, (int)j});
        
            }
        }   
    }

}

template <class filament_type>
vector<filament_type *>* filament_ensemble<filament_type>::get_network()
{
    return &network;
}

//given motor head position, return a map between  
//  the INDICES (i.e., {i, j} where i is the filament index and j is the link index)
//  and their corresponding DISTANCES to the link at that distance 
//NOTE: currently only looks for Link's IN the motor head's quadrant (determined by floor[pos/(fov)*nq])
template <class filament_type>
map<array<int,2>,double> filament_ensemble<filament_type>::get_dist(double x, double y)
{
    array<int, 2> motor_quad = {int(floor(x/fov[0]*nq[0])), int(floor(y/fov[1]*nq[1]))};
    map<array<int, 2>, double> t_map;
    if(!quad_fils[motor_quad].empty())
    {
        for (unsigned int j = 0; j < quad_fils[motor_quad].size(); j++)
            t_map[quad_fils[motor_quad][j]] = network[quad_fils[motor_quad][j][0]]->get_link(quad_fils[motor_quad][j][1])->get_distance(x,y);
    
    }
    return t_map;
}

template <class filament_type>
array<double,2> filament_ensemble<filament_type>::get_intpoints(int fil, int link, double xp, double yp)
{
    return network[fil]->get_link(link)->get_intpoint(xp,yp);
}

template <class filament_type> 
double filament_ensemble<filament_type>::get_xcm(int fil, int link)
{
    return network[fil]->get_link(link)->get_xcm();
}

template <class filament_type> 
double filament_ensemble<filament_type>::get_ycm(int fil, int link)
{
    return network[fil]->get_link(link)->get_ycm();
}

template <class filament_type> 
double filament_ensemble<filament_type>::get_angle(int fil, int link)
{
    return network[fil]->get_link(link)->get_angle();
}

template <class filament_type> 
double filament_ensemble<filament_type>::get_llength(int fil, int link)
{
    return network[fil]->get_link(link)->get_length();
}

template <class filament_type>
array<double,2> filament_ensemble<filament_type>::get_start(int fil, int link)
{
    return {network[fil]->get_link(link)->get_hx()[0] , network[fil]->get_link(link)->get_hy()[0]};
}

template <class filament_type>
array<double,2> filament_ensemble<filament_type>::get_end(int fil, int link)
{
    return {network[fil]->get_link(link)->get_hx()[1] , network[fil]->get_link(link)->get_hy()[1]};
}

template <class filament_type>
array<double,2> filament_ensemble<filament_type>::get_force(int fil, int actin)
{
    return network[fil]->get_actin(actin)->get_force();
}

template <class filament_type>
array<double,2> filament_ensemble<filament_type>::get_direction(int fil, int link)
{
    return network[fil]->get_link(link)->get_direction();
}

template <class filament_type> 
void filament_ensemble<filament_type>::set_straight_filaments(bool is_straight)
{
    straight_filaments = is_straight;
}

template <class filament_type> 
void filament_ensemble<filament_type>::update_positions(double t)
{
    for (unsigned int f = 0; f < network.size(); f++)
    {
        network[f]->update_positions(t);
    }

}

template <class filament_type> 
void filament_ensemble<filament_type>::write_actins(ofstream& fout)
{
    for (unsigned int i=0; i<network.size(); i++) {
        fout<<network[i]->write_actins();
    } 
}

template <class filament_type> 
void filament_ensemble<filament_type>::write_links(ofstream& fout)
{
    for (unsigned int i=0; i<network.size(); i++) {
        fout<<network[i]->write_links();
    } 
}

template <class filament_type> 
void filament_ensemble<filament_type>::write_thermo(ofstream& fout){
    for (unsigned int f = 0; f < network.size(); f++)
        fout<<network[f]->write_thermo();
    
}

template <class filament_type> 
void filament_ensemble<filament_type>::set_shear_rate(double g)
{
    for (unsigned int f = 0; f < network.size(); f++)
    {
        gamma = g;
    }
}


template <class filament_type> 
void filament_ensemble<filament_type>::update_shear(){
    
    for (unsigned int f = 0; f < network.size(); f++)
    {
        network[f]->update_shear();
    }
}

template <class filament_type> 
void filament_ensemble<filament_type>::print_filament_thermo(){
    
    for (unsigned int f = 0; f < network.size(); f++)
    {
        cout<<"\nF"<<f<<"\t:";
        network[f]->print_thermo();
    }

}

template <class filament_type> 
void filament_ensemble<filament_type>::print_network_thermo(){
    double KE=0, PE=0, TE=0;
    for (unsigned int f = 0; f < network.size(); f++)
    {
        KE += network[f]->get_kinetic_energy();
        PE += network[f]->get_potential_energy();
        TE += network[f]->get_total_energy();
    }
    cout<<"\nAll Fs\t:\tKE = "<<KE<<"\tPE = "<<PE<<"\tTE = "<<TE;
}


template <class filament_type> 
bool filament_ensemble<filament_type>::is_polymer_start(int fil, int actin){

    return !(actin);

}

template <class filament_type> 
void filament_ensemble<filament_type>::set_ld(double length){
    ld = length;
}

template <class filament_type> 
void filament_ensemble<filament_type>::set_fov(double fovx, double fovy){
    fov[0] = fovx;
    fov[1] = fovy;
}

template <class filament_type> 
void filament_ensemble<filament_type>::set_nq(double nqx, double nqy){
    nq[0] = nqx;
    nq[1] = nqy;
}

template <class filament_type> 
void filament_ensemble<filament_type>::set_visc(double nu){
    visc = nu;
}

template <class filament_type> 
void filament_ensemble<filament_type>::update_forces(int f_index, int a_index, double f1, double f2)
{
    network[f_index]->update_forces(a_index, f1,f2);
}

template <class filament_type> 
vector<int> filament_ensemble<filament_type>::get_broken(){
    return broken_filaments;
}

template <class filament_type> 
void filament_ensemble<filament_type>::clear_broken(){
    broken_filaments.clear();
}

// Update bending forces between monomers
template <class filament_type>
void filament_ensemble<filament_type>::update_bending(){
    
    for (unsigned int f = 0; f < network.size(); f++)
    {
        network[f]->update_bending();
    }
}

template <class filament_type>
void filament_ensemble<filament_type>::update_stretching(){
    
    vector<filament_type *> newfilaments;
    int s = network.size(); //keep it to one fracture per filament per timestep, or things get messy
    for (int f = 0; f < s; f++)
    {
        newfilaments = network[f]->update_stretching();
        
        if (newfilaments.size() > 0){ //fracture event occured

            cout<<"\n\tDEBUG: fracturing filament : "<<f;
            filament * broken = network[f];     //store a pointer to the broken filament to delete it with
            network[f] = newfilaments[0];       //replace that pointer with one of the new filaments
            
            if (newfilaments.size() == 2) network.push_back(newfilaments[1]); //add the second filament to the top of the stack
        
            broken_filaments.push_back(f);      // record the index, for automatic motor detachment
            delete broken;                      // delete the old filament
            
        }

    }
}

template <class filament_type>
void filament_ensemble<filament_type>::update(double t){
    
    this->update_shear();
    this->update_stretching();
    this->update_bending();
    this->update_positions(t);
    this->quad_update();

}

////////////////////////////////////////
///SPECIFIC FILAMENT IMPLEMENTATIONS////
////////////////////////////////////////

ATfilament_ensemble::ATfilament_ensemble(double density, array<double,2> myfov, array<int,2> mynq, double delta_t, double temp,
        double rad, double vis, int nactins, double link_len, vector<double *> pos_sets, double stretching, double bending, 
        double frac_force, string bc, double seed) {
    
    fov = myfov;
    nq = mynq;

    view[0] = 1;//(fov[0] - 2*nactins*link_len)/fov[0];
    view[1] = 1;//(fov[1] - 2*nactins*link_len)/fov[1];

    rho=density;
    visc=vis;
    ld=rad;//rng_n(len,1.0);
    link_ld = link_len;
    npolymer=int(ceil(density*fov[0]*fov[1]) / nactins);
    dt = delta_t;
    temperature = temp;

    if (seed == -1){
        straight_filaments = true;
    }else{
        srand(seed);
    }

    cout<<"DEBUG: Number of filament:"<<npolymer<<"\n";
    cout<<"DEBUG: Number of monomers per filament:"<<nactins<<"\n"; 
    cout<<"DEBUG: Monomer Length:"<<ld<<"\n"; 
    
    int s = pos_sets.size();
    double x0, y0, phi0;
    for (int i=0; i<npolymer; i++) {
        if ( i < s ){
            network.push_back(new filament({pos_sets[i][0], pos_sets[i][1], pos_sets[i][2]}, nactins, fov, nq,
                        visc, dt, temp, straight_filaments, ld, link_ld, stretching, bending, frac_force, bc) );
        }else{
            x0 = rng(-0.5*(view[0]*fov[0]),0.5*(view[0]*fov[0])); 
            y0 = rng(-0.5*(view[1]*fov[1]),0.5*(view[1]*fov[1]));
            phi0 =  rng(0, 2*pi);
            network.push_back(new filament({x0,y0,phi0}, nactins, fov, nq, visc, dt, temp, straight_filaments, ld, link_ld, stretching, bending, frac_force, bc) );
        }
    }
}

baoab_filament_ensemble::baoab_filament_ensemble(double density, array<double,2> myfov, array<int,2> mynq, double delta_t, double temp,
        double rad, double vis, int nactins, double link_len, vector<double *> pos_sets, double stretching, double bending, 
        double frac_force, string bc, double seed) {
    
    fov = myfov;
    nq = mynq;

    view[0] = 1;//(fov[0] - 2*nactins*link_len)/fov[0];
    view[1] = 1;//(fov[1] - 2*nactins*link_len)/fov[1];

    rho=density;
    visc=vis;
    ld=rad;//rng_n(len,1.0);
    link_ld = link_len;
    npolymer=int(ceil(density*fov[0]*fov[1]) / nactins);
    dt = delta_t;
    temperature = temp;

    if (seed == -1){
        straight_filaments = true;
    }else{
        srand(seed);
    }

    cout<<"DEBUG: Number of filament:"<<npolymer<<"\n";
    cout<<"DEBUG: Number of monomers per filament:"<<nactins<<"\n"; 
    cout<<"DEBUG: Monomer Length:"<<ld<<"\n"; 
    
    int s = pos_sets.size();
    double x0, y0, phi0;
    for (int i=0; i<npolymer; i++) {
        if ( i < s ){
            network.push_back(new baoab_filament({pos_sets[i][0], pos_sets[i][1], pos_sets[i][2]}, nactins, fov, nq,
                        visc, dt, temp, straight_filaments, ld, link_ld, stretching, bending, frac_force, bc) );
        }else{
            x0 = rng(-0.5*(view[0]*fov[0]),0.5*(view[0]*fov[0])); 
            y0 = rng(-0.5*(view[1]*fov[1]),0.5*(view[1]*fov[1]));
            phi0 =  rng(0, 2*pi);
            network.push_back(new baoab_filament({x0,y0,phi0}, nactins, fov, nq, visc, dt, temp, straight_filaments, 
                        ld, link_ld, stretching, bending, frac_force, bc) );
        }
    }
}

void baoab_filament_ensemble::update_velocities_B()
{
    for (unsigned int f = 0; f < network.size(); f++) 
        network[f]->update_velocities_B();
}
    
void baoab_filament_ensemble::update_velocities_O(double t)
{
    for (unsigned int f = 0; f < network.size(); f++) 
        network[f]->update_velocities_O(t);
}

void baoab_filament_ensemble::update(double t)
{

    this->update_velocities_B();
    this->update_positions(t);

    this->update_shear();
    this->update_stretching();
    this->update_bending();
    
    this->update_velocities_O(t);
    this->update_positions(t);
    this->update_velocities_B();
    
    this->quad_update();

}

lammps_filament_ensemble::lammps_filament_ensemble(double density, array<double,2> myfov, array<int,2> mynq, double delta_t, double temp,
        double rad, double vis, int nactins, double link_len, vector<double *> pos_sets, double stretching, double bending, 
        double frac_force, string bc, double seed) {
    
    fov = myfov;
    nq = mynq;

    view[0] = 1;//(fov[0] - 2*nactins*link_len)/fov[0];
    view[1] = 1;//(fov[1] - 2*nactins*link_len)/fov[1];

    rho=density;
    visc=vis;
    ld=rad;//rng_n(len,1.0);
    link_ld = link_len;
    npolymer=int(ceil(density*fov[0]*fov[1]) / nactins);
    dt = delta_t;
    temperature = temp;

    if (seed == -1){
        straight_filaments = true;
    }else{
        srand(seed);
    }

    cout<<"DEBUG: Number of filament:"<<npolymer<<"\n";
    cout<<"DEBUG: Number of monomers per filament:"<<nactins<<"\n"; 
    cout<<"DEBUG: Monomer Length:"<<ld<<"\n"; 
    
    int s = pos_sets.size();
    double x0, y0, phi0;
    for (int i=0; i<npolymer; i++) {
        if ( i < s ){
            network.push_back(new lammps_filament({pos_sets[i][0], pos_sets[i][1], pos_sets[i][2]}, nactins, fov, nq,
                        visc, dt, temp, straight_filaments, ld, link_ld, stretching, bending, frac_force, bc) );
        }else{
            x0 = rng(-0.5*(view[0]*fov[0]),0.5*(view[0]*fov[0])); 
            y0 = rng(-0.5*(view[1]*fov[1]),0.5*(view[1]*fov[1]));
            phi0 =  rng(0, 2*pi);
            network.push_back(new lammps_filament({x0,y0,phi0}, nactins, fov, nq, visc, dt, temp, straight_filaments, 
                        ld, link_ld, stretching, bending, frac_force, bc) );
        }
    }
    double mass_density = 2.6e-14; //miligram / micron
    this->set_mass(mass_density*link_len); 
}

void lammps_filament_ensemble::set_mass(double m)
{
    for (unsigned int f = 0; f < network.size(); f++) 
        network[f]->set_mass(m);
}

void lammps_filament_ensemble::update_brownian()
{
    for (unsigned int f = 0; f < network.size(); f++) 
        network[f]->update_brownian();
}

void lammps_filament_ensemble::update_drag()
{
    for (unsigned int f = 0; f < network.size(); f++) 
        network[f]->update_drag();
}

void lammps_filament_ensemble::update(double t){
    
    this->update_shear();
    this->update_stretching();
    this->update_bending();
    
    if (t > dt*100000) this->update_brownian();
    this->update_drag();
    
    this->update_positions(t);
    this->quad_update();

}

template class filament_ensemble<filament>;
template class filament_ensemble<baoab_filament>;
template class filament_ensemble<lammps_filament>;
