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

filament_ensemble::filament_ensemble(){}

filament_ensemble::filament_ensemble(double density, double fovx, double fovy, int nx, int ny, double delta_t, double temp,
        double len, double vis, int nrods, double link_len, std::vector<double *> pos_sets, double stretching, double bending, 
        double frac_force, double seed) {
    
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

    std::cout<<"DEBUG: Number of filament:"<<npolymer<<"\n";
    std::cout<<"DEBUG: Number of monomers per filament:"<<nrods<<"\n"; 
    std::cout<<"DEBUG: Monomer Length:"<<ld<<"\n"; 
    
    for (int i=0; i<npolymer; i++) {
        int s = pos_sets.size();
        if ( i < s){
            network.push_back(new filament(pos_sets[i][0], pos_sets[i][1], pos_sets[i][2], nrods, fov[0], fov[1], nq[0], nq[1],
                        visc, dt, temp, straight_filaments, ld, link_ld, stretching, bending, frac_force) );
        }else{
            network.push_back(new filament(rng(-0.5*(view[0]*fov[0]),0.5*(view[0]*fov[0])), rng(-0.5*(view[1]*fov[1]),0.5*(view[1]*fov[1])), rng(0, 2*pi),
                        nrods, fov[0], fov[1], nq[0], nq[1],
                        visc, dt, temp, straight_filaments, ld, link_ld, stretching, bending, frac_force) );
        }
    }
}

filament_ensemble::~filament_ensemble(){ 
    std::cout<<"DELETING filament_ensemble\n";
    
    int s = network.size();
    for (int i = 0; i < s; i++){
        delete network[i];
    }
    
    network.clear();
};

void filament_ensemble::quad_update()
{
    quad_fils.clear();
    std::vector<std::vector<std::vector<int> > > all_tmp_quads;
    std::vector<std::vector<int> > tmp_quads;
    std::vector<int> index;

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

std::vector<filament *>* filament_ensemble::get_network()
{
    return &network;
}

//given motor quadrant, return the indices and the distances to all actin filaments in the neighboring quadrants 
std::map<std::vector<int>,double> filament_ensemble::get_dist(double x, double y)
{
    int nxx=int(floor(x/fov[0]*nq[0]));
    int nyy=int(floor(y/fov[1]*nq[1]));
    std::vector<int> index;
    t_map.clear();
    if(!quad_fils[nxx][nyy].empty())
    {
        for(std::vector<std::vector<int> >::iterator it=quad_fils[nxx][nyy].begin(); it<quad_fils[nxx][nyy].end(); it++)
        {   
            index = *it;
            t_map[*it] = network[ index[0] ]->get_rod( index[1])->get_distance(x,y);
        }
    }
    return t_map;
}

double* filament_ensemble::get_intpoints(int fil, int rod, double xp, double yp)
{
    return network[fil]->get_rod(rod)->get_intpoint(xp,yp);
}

double filament_ensemble::get_xcm(int fil, int rod)
{
    return network[fil]->get_rod(rod)->get_xcm();
}

double filament_ensemble::get_ycm(int fil, int rod)
{
    return network[fil]->get_rod(rod)->get_ycm();
}

double filament_ensemble::get_angle(int fil, int rod)
{
    return network[fil]->get_rod(rod)->get_angle();
}

double filament_ensemble::get_alength(int fil, int rod)
{
    return network[fil]->get_rod(rod)->get_length();
}

double * filament_ensemble::get_start(int fil, int rod)
{
    return network[fil]->get_rod(rod)->get_start();
}

double * filament_ensemble::get_end(int fil, int rod)
{
    return network[fil]->get_rod(rod)->get_end();
}

double * filament_ensemble::get_forces(int fil, int rod)
{
    return network[fil]->get_rod(rod)->get_forces();
}

void filament_ensemble::set_straight_filaments(bool is_straight)
{
    straight_filaments = is_straight;
}

void filament_ensemble::update(double t)
{
    for (unsigned int f = 0; f < network.size(); f++)
    {
        network[f]->update(t);
    }

}

void filament_ensemble::write(std::ofstream& fout)
{
    for (unsigned int i=0; i<network.size(); i++) {
        fout<<network.at(i)->write();
    } 
}


// Update bending forces between monomers
void filament_ensemble::update_bending(){
    
    for (unsigned int f = 0; f < network.size(); f++)
    {
        network[f]->update_bending();
    }
}

void filament_ensemble::update_stretching(){

    std::vector<filament *> newfilaments;
    for (unsigned int f = 0; f < network.size(); f++)
    {
        newfilaments = network[f]->update_stretching();
        
        if (newfilaments){ //fracture event occured
            network.erase(network.begin() + f);
            network.push_back(newfilaments[0]);
            network.push_back(newfilaments[1]);
        }

    }
}

void filament_ensemble::set_shear_rate(double g)
{
    for (unsigned int f = 0; f < network.size(); f++)
    {
        gamma = g;
    }
}


void filament_ensemble::update_shear(){
    
    for (unsigned int f = 0; f < network.size(); f++)
    {
        network[f]->update_shear();
    }
}

bool filament_ensemble::is_polymer_start(int fil, int rod){

    return !(rod);

}

void filament_ensemble::set_ld(double length){
    ld = length;
}

void filament_ensemble::set_fov(double fovx, double fovy){
    fov[0] = fovx;
    fov[1] = fovy;
}

void filament_ensemble::set_nq(double nqx, double nqy){
    nq[0] = nqx;
    nq[1] = nqy;
}

void filament_ensemble::set_visc(double nu){
    visc = nu;
}

