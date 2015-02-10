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
    cout<<"DELETING filament_ensemble\n";
    
    int s = network.size();
    for (int i = 0; i < s; i++){
        delete network[i];
    }
    
    network.clear();
};

template <class filament_type> 
void filament_ensemble<filament_type>::quad_update()
{
    quad_fils.clear();
    vector<vector<vector<int> > > all_tmp_quads;
    vector<vector<int> > tmp_quads;
    vector<int> index;

    for (unsigned int i=0; i<network.size(); i++) {
        
        all_tmp_quads = network[i]->get_quadrants(); //quadrants for every rod on the filament
        
        for (unsigned int j=0; j<all_tmp_quads.size(); j++){
            
            tmp_quads = all_tmp_quads[j]; //quadrants of the jth rod on the filament
            index.push_back(i);
            index.push_back(j);

            for (unsigned int xindex=0; xindex<tmp_quads[0].size(); xindex++) {
                for (unsigned int yindex=0; yindex<tmp_quads[1].size(); yindex++) {
                    quad_fils[tmp_quads[0][xindex]][tmp_quads[1][yindex]].push_back(index);
                }
            }

            index.clear();

        }
    }
}

template <class filament_type>
vector<filament_type *>* filament_ensemble<filament_type>::get_network()
{
    return &network;
}

//given motor quadrant, return the indices and the distances to all actin filaments in the neighboring quadrants 
template <class filament_type>
map<vector<int>,double> filament_ensemble<filament_type>::get_dist(double x, double y)
{
    int nxx=int(floor(x/fov[0]*nq[0]));
    int nyy=int(floor(y/fov[1]*nq[1]));
    vector<int> index;
    t_map.clear();
    if(!quad_fils[nxx][nyy].empty())
    {
        for(vector<vector<int> >::iterator it=quad_fils[nxx][nyy].begin(); it<quad_fils[nxx][nyy].end(); it++)
        {   
            index = *it;
            t_map[*it] = network[ index[0] ]->get_rod( index[1])->get_distance(x,y);
        }
    }
    return t_map;
}

template <class filament_type>
double* filament_ensemble<filament_type>::get_intpoints(int fil, int rod, double xp, double yp)
{
    return network[fil]->get_rod(rod)->get_intpoint(xp,yp);
}

template <class filament_type> 
double filament_ensemble<filament_type>::get_xcm(int fil, int rod)
{
    return network[fil]->get_rod(rod)->get_xcm();
}

template <class filament_type> 
double filament_ensemble<filament_type>::get_ycm(int fil, int rod)
{
    return network[fil]->get_rod(rod)->get_ycm();
}

template <class filament_type> 
double filament_ensemble<filament_type>::get_angle(int fil, int rod)
{
    return network[fil]->get_rod(rod)->get_angle();
}

template <class filament_type> 
double filament_ensemble<filament_type>::get_alength(int fil, int rod)
{
    return network[fil]->get_rod(rod)->get_length();
}

template <class filament_type>
double * filament_ensemble<filament_type>::get_start(int fil, int rod)
{
    return network[fil]->get_rod(rod)->get_start();
}

template <class filament_type>
double * filament_ensemble<filament_type>::get_end(int fil, int rod)
{
    return network[fil]->get_rod(rod)->get_end();
}

template <class filament_type>
double * filament_ensemble<filament_type>::get_forces(int fil, int rod)
{
    return network[fil]->get_rod(rod)->get_forces();
}

template <class filament_type>
double * filament_ensemble<filament_type>::get_direction(int fil, int rod)
{
    return network[fil]->get_rod(rod)->get_direction();
}

template <class filament_type> 
void filament_ensemble<filament_type>::set_straight_filaments(bool is_straight)
{
    straight_filaments = is_straight;
}

template <class filament_type> 
void filament_ensemble<filament_type>::update(double t)
{
    for (unsigned int f = 0; f < network.size(); f++)
    {
        network[f]->update(t);
    }

}

template <class filament_type> 
void filament_ensemble<filament_type>::write_rods(ofstream& fout)
{
    for (unsigned int i=0; i<network.size(); i++) {
        fout<<network[i]->write_rods();
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
bool filament_ensemble<filament_type>::is_polymer_start(int fil, int rod){

    return !(rod);

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
void filament_ensemble<filament_type>::update_forces(int f_index, int r_index, double f1, double f2, double f3)
{
    network[f_index]->update_forces(r_index, f1,f2,f3);
}

////////////////////////////////////////
///SPECIFIC FILAMENT IMPLEMENTATIONS////
////////////////////////////////////////

ATfilament_ensemble::ATfilament_ensemble(double density, double fovx, double fovy, int nx, int ny, double delta_t, double temp,
        double len, double vis, int nrods, double link_len, vector<double *> pos_sets, double stretching, double bending, 
        double frac_force, string bc, double seed) {
    
    view[0] = (fovx - 2*nrods*len)/fovx;
    view[1] = (fovy - 2*nrods*len)/fovy;

    fov[0]=fovx;
    fov[1]=fovy;
    nq[0]=nx;
    nq[1]=ny;
    rho=density;
    visc=vis;
    ld=len;//rng_n(len,1.0);
    link_ld = link_len;
    npolymer=int(ceil(density*fov[0]*fov[1]) / nrods);
    dt = delta_t;
    temperature = temp;

    if (seed == -1){
        straight_filaments = true;
    }else{
        srand(seed);
    }

    cout<<"DEBUG: Number of filament:"<<npolymer<<"\n";
    cout<<"DEBUG: Number of monomers per filament:"<<nrods<<"\n"; 
    cout<<"DEBUG: Monomer Length:"<<ld<<"\n"; 
    
    int s = pos_sets.size();
    
    for (int i=0; i<npolymer; i++) {
        if ( i < s ){
            network.push_back(new filament(pos_sets[i][0], pos_sets[i][1], pos_sets[i][2], nrods, fov[0], fov[1], nq[0], nq[1],
                        visc, dt, temp, straight_filaments, ld, link_ld, stretching, bending, frac_force, bc) );
        }else{
            network.push_back(new filament(rng(-0.5*(view[0]*fov[0]),0.5*(view[0]*fov[0])), 
                        rng(-0.5*(view[1]*fov[1]),0.5*(view[1]*fov[1])), rng(0, 2*pi),
                        nrods, fov[0], fov[1], nq[0], nq[1],
                        visc, dt, temp, straight_filaments, ld, link_ld, stretching, bending, frac_force, bc) );
        }
    }
}

// Update bending forces between monomers
void ATfilament_ensemble::update_bending(){
    
    for (unsigned int f = 0; f < network.size(); f++)
    {
        network[f]->update_bending();
    }
}

void ATfilament_ensemble::update_stretching(){
    
    vector<filament *> newfilaments;
    for (unsigned int f = 0; f < network.size(); f++)
    {
        newfilaments = network[f]->update_stretching();
        
        if (newfilaments.size() > 0){ //fracture event occured
            filament * broken = network[f];
            network.erase(network.begin() + f);
            network.push_back(newfilaments[0]);
            network.push_back(newfilaments[1]);
            delete broken;
            break; //avoid infinite loops
        }

    }
}

DLfilament_ensemble::DLfilament_ensemble(double density, double fovx, double fovy, int nx, int ny, double delta_t, double temp,
        double len, double vis, int nrods, double link_len, vector<double *> pos_sets, double stretching, double bending, 
        double frac_force, double bending_frac_force, string bc, double seed) {
    
    view[0] = (fovx - 2*nrods*len)/fovx;
    view[1] = (fovy - 2*nrods*len)/fovy;

    fov[0]=fovx;
    fov[1]=fovy;
    nq[0]=nx;
    nq[1]=ny;
    rho=density;
    visc=vis;
    ld=len;//rng_n(len,1.0);
    link_ld = link_len;
    npolymer=int(ceil(density*fov[0]*fov[1]) / nrods);
    dt = delta_t;
    temperature = temp;

    if (seed == -1){
        straight_filaments = true;
    }else{
        srand(seed);
    }

    int s = pos_sets.size();
    
    for (int i=0; i<npolymer; i++) {
        if ( i < s ){
            network.push_back(new DLfilament(pos_sets[i][0], pos_sets[i][1], pos_sets[i][2], nrods, fov[0], fov[1], nq[0], nq[1],
                        visc, dt, temp, straight_filaments, ld, link_ld, stretching, bending, frac_force, bending_frac_force, bc) );
        }else{
            network.push_back(new DLfilament(rng(-0.5*(view[0]*fov[0]),0.5*(view[0]*fov[0])), 
                        rng(-0.5*(view[1]*fov[1]),0.5*(view[1]*fov[1])), rng(0, 2*pi),
                        nrods, fov[0], fov[1], nq[0], nq[1],
                        visc, dt, temp, straight_filaments, ld, link_ld, stretching, bending, frac_force, bending_frac_force, bc) );
        }
    }
}

void DLfilament_ensemble::update_stretching(){

    vector<DLfilament *> newfilaments;
    for (unsigned int f = 0; f < network.size(); f++)
    {
        newfilaments = network[f]->update_stretching();
        
        if (newfilaments.size() > 0){ //fracture event occured
            DLfilament * broken = network[f];
            network.erase(network.begin() + f);
            network.push_back(newfilaments[0]);
            network.push_back(newfilaments[1]);
            delete broken; 
            break; //avoid infinite loops
        }

    }
}

void DLfilament_ensemble::update_bending(){

    vector<DLfilament *> newfilaments;
    for (unsigned int f = 0; f < network.size(); f++)
    {
        newfilaments = network[f]->update_bending();
        
        if (newfilaments.size() > 0){ //fracture event occured
            network.erase(network.begin() + f);
            network.push_back(newfilaments[0]);
            network.push_back(newfilaments[1]);
            break; //avoid infinite loops
        }

    }
}

void DLfilament_ensemble::set_bending_linear(){
    for (unsigned int f = 0; f < network.size(); f++)
        network[f]->set_bending_linear();
}

template class filament_ensemble<filament>;
template class filament_ensemble<DLfilament>;
