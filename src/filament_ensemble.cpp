/*-------------------------------------------------------------------
 filament_ensemble.cpp : container class for filaments
 
 Copyright (C) 2016 
 Created by: Simon Freedman, Shiladitya Banerjee, Glen Hocky, Aaron Dinner
 Contact: dinner@uchicago.edu

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version. See ../LICENSE for details. 
-------------------------------------------------------------------*/

#include "globals.h"
#include "filament_ensemble.h"
//bead network class

filament_ensemble::filament_ensemble(){}

 
filament_ensemble::~filament_ensemble(){ 
    cout<<"DELETING FILAMENT_ENSEMBLE\n";
    
    int s = network.size();
    
    for (int x = 0; x < nq[0]; x++){
        for (int y = 0; y < nq[1]; y++){
            delete springs_per_quad[x]->at(y);
        }
        delete springs_per_quad[x];
        //delete n_springs_per_quad[x];
    }
    
    for (int i = 0; i < s; i++){
        delete network[i];
    }
    
}


vector<filament *>* filament_ensemble::get_network()
{
    return &network;
}


filament * filament_ensemble::get_filament(int index)
{
    return network[index];
}


void filament_ensemble::turn_quads_off()
{
    quad_off_flag = true;
}


void filament_ensemble::nlist_init_serial()
{
    for (int x = 0; x < nq[0]; x++){
        springs_per_quad.push_back(new vector< vector<array<int, 2> >* >(nq[1]));   
        for (int y = 0; y < nq[1]; y++){
            springs_per_quad[x]->at(y) = new vector<array<int, 2> >();
        }
    }
}

 
void filament_ensemble::quad_update_serial()
{
    int n_quads, net_sz = int(network.size());
    vector<vector<array<int, 2> > > q;
    int x, y;

    //initialize all quadrants to have no springs
    for (x = 0; x < nq[0]; x++){
        for (y = 0; y < nq[1]; y++){
            springs_per_quad[x]->at(y)->clear();
        }
    }
    
    for (int f = 0; f < net_sz; f++){
        q = network[f]->get_quadrants();
        for (int l = 0; l < network[f]->get_nsprings(); l++){
            n_quads = int(q[l].size());
            for (int i = 0; i < n_quads; i++){
                x = q[l][i][0];
                y = q[l][i][1];
                springs_per_quad[x]->at(y)->push_back({{f,l}});
                
            }
        }
    }
}

//given a motor position, and a quadrant
//update the map of {{f, l}} -- > dist
void filament_ensemble::update_dist_map(set<pair<double, array<int,2>>>& t_map, const array<int, 2>& mq, double x, double y){
    
    array<int, 2> fl;
    double dist_sq;
    
    for (int i = 0; i < int(springs_per_quad[mq[0]]->at(mq[1])->size()); i++){

        fl = springs_per_quad[mq[0]]->at(mq[1])->at(i); //fl  = {{filament_index, spring_index}}

        if (fls.find(fl) == fls.end()){
            network[fl[0]]->get_spring(fl[1])->calc_intpoint(network[fl[0]]->get_BC(), delrx, x, y); //calculate the point on the spring closest to (x,y)
            dist_sq = network[fl[0]]->get_spring(fl[1])->get_distance_sq(network[fl[0]]->get_BC(), delrx, x, y); //store the distance to that point
            //cout<<"\nDEBUG : dist = "<<dist;

            t_map.insert(pair<double, array<int, 2> >(dist_sq, fl));
            fls.insert(fl);
        }
    }

}

//given motor head position, return a map between  
//  the INDICES (i.e., {i, j} for the j'th spring of the i'th filament)
//  and their corresponding DISTANCES to the spring at that distance 

set<pair<double, array<int, 2>>> filament_ensemble::get_dist(double x, double y)
{
    fls.clear();
    set<pair<double, array<int, 2>>> t_map;
    
    int mqx = coord2quad(fov[0], nq[0], x);
    int mqy = coord2quad(fov[1], nq[1], y);
    
    update_dist_map(t_map, {{mqx, mqy}}, x, y);
    return t_map;
}


set<pair<double, array<int,2>>> filament_ensemble::get_dist_all(double x, double y)
{
    set<pair<double, array<int,2>>> t_map;
    double dist_sq=0;
    for (int f = 0; f < int(network.size()); f++){
        for (int l=0; l < network[f]->get_nsprings(); l++){
                network[f]->get_spring(l)->calc_intpoint(network[f]->get_BC(), delrx, x, y); //calculate the point on the spring closest to (x,y)
                dist_sq = network[f]->get_spring(l)->get_distance_sq(network[f]->get_BC(), delrx, x, y); //store the distance to that point
                // t_map[dist] = {f,l}; 
                t_map.insert(pair<double, array<int, 2>>(dist_sq, {{f, l}}));
        }
    }
    
    return t_map;
}

double filament_ensemble::get_llength(int fil, int spring)
{
    return network[fil]->get_spring(spring)->get_length();
}


array<double,2> filament_ensemble::get_start(int fil, int spring)
{
    return {{network[fil]->get_spring(spring)->get_hx()[0] , network[fil]->get_spring(spring)->get_hy()[0]}};
}


array<double,2> filament_ensemble::get_end(int fil, int spring)
{
    return {{network[fil]->get_spring(spring)->get_hx()[1] , network[fil]->get_spring(spring)->get_hy()[1]}};
}


array<double,2> filament_ensemble::get_force(int fil, int bead)
{
    return network[fil]->get_bead(bead)->get_force();
}


array<double,2> filament_ensemble::get_direction(int fil, int spring)
{
    return network[fil]->get_spring(spring)->get_direction();
}

 
void filament_ensemble::set_straight_filaments(bool is_straight)
{
    straight_filaments = is_straight;
}

 
void filament_ensemble::update_positions()
{
    int net_sz = int(network.size());
    for (int f = 0; f < net_sz; f++)
    {
        //if (f==0) cout<<"\nDEBUG: update_positions: using "<<omp_get_num_threads()<<" cores";  
        network[f]->update_positions();
    }

}

 
void filament_ensemble::update_positions_range(int lo, int hi)
{
    for (unsigned int f = 0; f < network.size(); f++)
    {
        network[f]->update_positions_range(lo, hi);
    }

}

 
void filament_ensemble::write_beads(ofstream& fout)
{
    for (unsigned int i=0; i<network.size(); i++) {
        fout<<network[i]->write_beads(i);
    } 
}

 
void filament_ensemble::write_springs(ofstream& fout)
{
    for (unsigned int i=0; i<network.size(); i++) {
        fout<<network[i]->write_springs(i);
    } 
}

 
void filament_ensemble::write_thermo(ofstream& fout){
    for (unsigned int f = 0; f < network.size(); f++)
        fout<<network[f]->write_thermo(f);
    
}

 
void filament_ensemble::set_shear_rate(double g)
{
    if (network.size() > 0)
        if (network[0]->get_nbeads() > 0)
            shear_speed = g*fov[1] / (2*network[0]->get_bead(0)->get_friction());

    for (unsigned int f = 0; f < network.size(); f++)
    {
        network[f]->set_shear(g);
    }
}

 
void filament_ensemble::set_y_thresh(double y)
{
    for (unsigned int f = 0; f < network.size(); f++) network[f]->set_y_thresh(y);
}

 
void filament_ensemble::update_delrx(double drx)
{
    //cout<<"\nDEBUG: SHEARING"; 
    delrx = drx;
    for (unsigned int f = 0; f < network.size(); f++)
    {
        network[f]->update_delrx(drx);
    }
}

 
void filament_ensemble::update_d_strain(double g)
{
    //cout<<"\nDEBUG: SHEARING"; 
    for (unsigned int f = 0; f < network.size(); f++)
    {
        network[f]->update_d_strain(g);
    }
}

 
void filament_ensemble::update_shear()
{
    //cout<<"\nDEBUG: SHEARING"; 
    for (unsigned int f = 0; f < network.size(); f++)
    {
        network[f]->update_shear(t);
    }
}

 
void filament_ensemble::print_filament_thermo(){
    
    for (unsigned int f = 0; f < network.size(); f++)
    {
        cout<<"\nF"<<f<<"\t:";
        network[f]->print_thermo();
    }

}

 
void filament_ensemble::update_energies(){
    pe_stretch = 0;
    pe_bend = 0;
    ke_vir = 0; 
    for (unsigned int f = 0; f < network.size(); f++)
    { 
        ke_vir += network[f]->get_kinetic_energy(); 
        pe_bend += network[f]->get_bending_energy();
        pe_stretch += network[f]->get_stretching_energy();  
    }
}

 
double filament_ensemble::get_stretching_energy(){
    return pe_stretch;
}

 
double filament_ensemble::get_bending_energy(){
    return pe_bend;
}

double filament_ensemble::get_kinetic_energy_vir(){ 
    return ke_vir; 
}

void filament_ensemble::print_network_thermo(){
    cout<<"\nAll Fs\t:\tKE = "<<ke_vir<<"\tPEs = "<<pe_stretch<<"\tPEb = "<<pe_bend<<"\tPEexv = "<<pe_exv<<"\tTE = "<<(ke_vir+pe_stretch+pe_bend+pe_exv);
}

 
void filament_ensemble::print_filament_lengths(){
    for (unsigned int f = 0; f < network.size(); f++)
    {
        cout<<"\nF"<<f<<" : "<<network[f]->get_end2end()<<" um";
    }
}
 
bool filament_ensemble::is_polymer_start(int fil, int bead){

    return !(bead);

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

 
void filament_ensemble::update_forces(int f_index, int a_index, double f1, double f2){
    network[f_index]->update_forces(a_index, f1,f2);
}

 
vector<int> filament_ensemble::get_broken(){
    return broken_filaments;
}

 
void filament_ensemble::clear_broken(){
    broken_filaments.clear();
}

 
int filament_ensemble::get_nbeads(){
    int tot = 0;
    for (unsigned int f = 0; f < network.size(); f++)
        tot += network[f]->get_nbeads();
    return tot;
}

 
int filament_ensemble::get_nsprings(){
    return this->get_nbeads() - network.size();
}

 
int filament_ensemble::get_nfilaments(){
    return network.size();
}

 
double filament_ensemble::get_delrx(){
    return delrx;
}

 
double filament_ensemble::get_bead_friction(){
    
    if (network.size() > 0)
        if (network[0]->get_nbeads() > 0)
            return network[0]->get_bead(0)->get_friction();
    
    return 0;
}

// Update bending forces between monomers

void filament_ensemble::update_bending()
{
    int net_sz = int(network.size());
    
    for (int f = 0; f < net_sz; f++)
    {
        //if (f==0) cout<<"\nDEBUG: update_bending: using "<<omp_get_num_threads()<<" cores";  
        network[f]->update_bending(t);
    }
}


void filament_ensemble::update_stretching(){
    
//    vector<filament *> newfilaments;
    int s = network.size(); //keep it to one fracture per filament per timestep, or things get messy
    for (int f = 0; f < s; f++)
    {
        //if (f==0) cout<<"\nDEBUG: update_stretching: using "<<omp_get_num_threads()<<" cores";  
        this->update_filament_stretching(f);
    }
}


void filament_ensemble::update_filament_stretching(int f){
    vector<filament *> newfilaments = network[f]->update_stretching(t);

    if (newfilaments.size() > 0){ //fracture event occured

        cout<<"\n\tDEBUG: fracturing filament : "<<f;
        filament * broken = network[f];     //store a pointer to the broken filament to delete it with
        network[f] = newfilaments[0];       //replace that pointer with one of the new filaments

        if (newfilaments.size() == 2) network.push_back(newfilaments[1]); //add the second filament to the top of the stack

        broken_filaments.push_back(f);      // record the index, for automatic motor detachment
        delete broken;                      // delete the old filament

    }
}


void filament_ensemble::set_shear_stop(double stopT){
    shear_stop = stopT; 
}


void filament_ensemble::set_shear_dt(double delT){
    shear_dt = delT; 
}


void filament_ensemble::update_int_forces()
{
    this->update_stretching();
    this->update_bending();
}

/* Overdamped Langevin Dynamics Integrator (Leimkuhler, 2013) */

void filament_ensemble::update()
{      
    int net_sz = network.size();
    // #pragma omp parallel for     

    pe_exv = 0; 

    if ( kexv > 0)
        this->update_spring_forces_from_quads();

    for (int f = 0; f < net_sz; f++){
        
        network[f]->update_length();
        this->update_filament_stretching(f);
        network[f]->update_bending(t);
        network[f]->update_positions();
    }
 
    this->update_energies(); 
    
    t += dt;
}

void filament_ensemble::update_spring_forces_from_quads()
{
    //This function loops through the quads of the systems and then loops through the filaments and springs described by the neighbor list in every quad.
    //Upon this looping, the pair interactions will be calculated according to the neighbor list. 
    //The loop will go through thee values of spring_per_quad() [x][y][i], where x: 0-nq[0], y: 0-nq[1], i: 0-m_springs_per_quad
    //The pairs found in he nieghbor list are then saved in a vector array. 
    //On subsequent loops, this value will be searched for in order to ensure no repeats in the force calculation.  
    
    array <int,2> spring_1; 
    array <int,2> spring_2; 
    int f1, f2, l1, l2, nsprings_at_quad; 
    double par1, par2;   
    int max_nsprings = network.size()*nsprings_per_fil_max;
    
    vector<vector<int>> int_lks (max_nsprings, vector<int> (max_nsprings, 0));    
    
    for(int x = 0; x < nq[0]; x++) 
    {   
        for(int y = 0; y < nq[1]; y++) 
        { 
            nsprings_at_quad = int(springs_per_quad[x]->at(y)->size());
            
            for(int i = 0; i < nsprings_at_quad; i++) 
            {   
                spring_1 = springs_per_quad[x]->at(y)->at(i); 

                for(int j = i+1; j < nsprings_at_quad; j++) 
                {   
                    spring_2 = springs_per_quad[x]->at(y)->at(j);

                    f1 = spring_1[0]; 
                    f2 = spring_2[0]; 
                    l1 = spring_1[1];
                    l2 = spring_2[1]; 

                    par1 = f1*(network[f1]->get_nsprings()) + l1; 
                    par2 = f2*(network[f2]->get_nsprings()) + l2;

                    if ( f1 == f2 && abs(l1 - l2) < 2 ) // adjacent springs would yield excluded volume interactions between the same bead
                        continue;
                    else
                    {
                        par1 = f1*nsprings_per_fil_max + l1; 
                        par2 = f2*nsprings_per_fil_max + l2;

                        if ( ! int_lks[par1][par2] )
                        {
                            int_lks[par1][par2] = 1; 
                            int_lks[par2][par1] = 1; 

                            this->update_force_between_filaments(f1, l1, f2, l2);
                        }
                    }
                }
            }
        }
    }   
    int_lks.clear(); 
}

void filament_ensemble::update_spring_forces(int f) 
{
    //This function loops through every filament and spring in the network and applies the force calulation under certain limits

    int net_sz = network.size();
    int lks_sz = network[f]->get_nsprings();
    int oth_lks_sz;

    for(int i = 0; i < lks_sz; i++) 
    {   
        for(int g = f+1; g < net_sz; g++)
        {
            oth_lks_sz = network[g]->get_nsprings(); 

            for(int j = 0; j < oth_lks_sz; j++) 
            {
                this->update_force_between_filaments(f, i, g, j); 
            }
        }
    } 

}
void filament_ensemble::update_force_between_filaments(double n1, double l1, double n2, double l2)
{ 
    //This function calculates the forces applied to the actin beads of a pair of filaments under certain limits. 
    //Here, we use distance of closest approach to describe the direction and magnitude of the forces. 

    array <double, 4> r_c; 
    array <double, 2> p1, p2, p3, p4; 
    array <double, 2> len, hx_1, hy_1, hx_2, hy_2, dist;  
    double b = (1/rmax); 
    double r=0, x1=0, y1=0, x2=0, y2=0, length=0, len1=0, len2=0, r_1=0, r_2=0, Fx1=0, Fy1=0, Fx2=0, Fy2=0; 
    int index; 
    bool intersect; 

    hx_1 = network[n1]->get_spring(l1)->get_hx(); 
    hy_1 = network[n1]->get_spring(l1)->get_hy(); 

    hx_2 = network[n2]->get_spring(l2)->get_hx();
    hy_2 = network[n2]->get_spring(l2)->get_hy();

    r_c[0] = network[n1]->get_spring(l1)->get_r_c(BC, delrx, hx_2[0], hy_2[0]);
    p1 = network[n1]->get_spring(l1)->get_point(); 

    r_c[1] = network[n1]->get_spring(l1)->get_r_c(BC, delrx, hx_2[1], hy_2[1]);
    p2 = network[n1]->get_spring(l1)->get_point(); 
    
    r_c[2] = network[n2]->get_spring(l2)->get_r_c(BC, delrx, hx_1[0], hy_1[0]);
    p3 = network[n2]->get_spring(l2)->get_point(); 

    r_c[3] = network[n2]->get_spring(l2)->get_r_c(BC, delrx, hx_1[1], hy_1[1]);
    p4 = network[n2]->get_spring(l2)->get_point(); 

    len[0] = network[n1]->get_spring(l1)->get_length(); 
    len[1] = network[n2]->get_spring(l2)->get_length(); 
    
    r = r_c[0];
    index = 0; 

    for(int k = 1; k < 4; k++){
        if(r_c[k] < r){
            r = r_c[k]; 
            index = k; 
        }
    }   

    spring *L2 = network[n2]->get_spring(l2);  
    intersect = network[n1]->get_spring(l1)->get_line_intersect(BC, delrx, L2); 

    if(r < rmax)
    {
        if( !intersect )
        { 
            if(index == 0)
            {  
                r = r_c[0]; 
                x1 = hx_2[0];
                y1 = hy_2[0];
                x2 = p1[0];
                y2 = p1[1];
                len1 = dist_bc(BC, (hx_1[0]-x2), (hy_1[0]-y2), fov[0], fov[1], delrx);
                length = len[0];         
            }
            else if(index == 1)
            {
                r = r_c[1]; 
                x1 = hx_2[1];
                y1 = hy_2[1];
                x2 = p2[0]; 
                y2 = p2[1];
                len1 = dist_bc(BC, (hx_1[0]-x2), (hy_1[0]-y2), fov[0], fov[1], delrx);
                length = len[0];
            }
            else if(index == 2)
            {
                r = r_c[2];   
                x1 = hx_1[0]; 
                y1 = hy_1[0]; 
                x2 = p3[0];
                y2 = p3[1];
                len1 = dist_bc(BC, (hx_2[0]-x2), (hy_2[0]-y2), fov[0], fov[1], delrx);
                length = len[1]; 
            }
            else if(index == 3)
            {
                r = r_c[3]; 
                x1 = hx_1[1];
                y1 = hy_1[1];
                x2 = p4[0];
                y2 = p4[1];
                len1 = dist_bc(BC, (hx_2[0]-x2), (hy_2[0]-y2), fov[0], fov[1], delrx);
                length = len[1]; 
            }

            dist = rij_bc(BC, (x2-x1), (y2-y1), fov[0], fov[1], delrx); 
            len2 = length - len1; 
            r_1 = (len2/length);
            r_2 = (len1/length);

            Fx1 = 2*kexv*dist[0]*b*((1/r) - b); 
            Fx2 = -Fx1;
            Fy1 = 2*kexv*dist[1]*b*((1/r) - b);
            Fy2 = -Fy1;

            pe_exv += kexv*pow((1-r*b),2);

            if(index == 0)
            {
                network[n1]->update_forces(l1, Fx1*r_1, Fy1*r_1);
                network[n1]->update_forces(l1+1, Fx1*r_2, Fy1*r_2);
                network[n2]->update_forces(l2, Fx2, Fy2);
            }
            else if(index == 1)
            {
                network[n1]->update_forces(l1, Fx1*r_1, Fy1*r_1);
                network[n1]->update_forces(l1+1, Fx1*r_2, Fy1*r_2);
                network[n2]->update_forces(l2+1, Fx2, Fy2);
            }
            else if(index == 2)
            {
                network[n2]->update_forces(l2, Fx1*r_1, Fy1*r_1);
                network[n2]->update_forces(l2+1, Fx1*r_2, Fy1*r_2);
                network[n1]->update_forces(l1, Fx2, Fy2);
            }
            else if(index == 3)
            {
                network[n2]->update_forces(l2, Fx1*r_1, Fy1*r_1);
                network[n2]->update_forces(l2+1, Fx1*r_2, Fy1*r_2);
                network[n1]->update_forces(l1+1, Fx2, Fy2);
            } 
        }
        else
        { 
            Fx1 = 2*kexv/(rmax*sqrt(2)); 
            Fx2 = -Fx1; 
            Fy1 = 2*kexv/(rmax*sqrt(2)); 
            Fy2 = -Fy1; 

            pe_exv += kexv*pow((1-r*b),2);   

            network[n1]->update_forces(l1, Fx1, Fy1); 
            network[n1]->update_forces(l1+1, Fx1, Fy1); 
            network[n2]->update_forces(l2, Fx2, Fy2); 
            network[n2]->update_forces(l2+1, Fx2, Fy2); 
        }   
    }
}

double filament_ensemble::get_exv_energy()
{
    return pe_exv; 
} 

void filament_ensemble::update_excluded_volume(int f)
{
    //For every filament bead on f, for every bead not on f, calculate the force between the two bead using the Jones potential, and update them ( maybe divide by half due to overcaluclations).	

    int net_sz = network.size();
    int act_sz = network[f]->get_nbeads();  
    //10^6 included to account for m to microm conversion
    double a = 0.004; 
    double b = 1/rmax; 
    double x1, x2, y1, y2, Fx1, Fx2, Fy1, Fy2, r, dx, dy; 


    for(int i = 0; i < act_sz; i++){
        for(int g = f+1; g < net_sz; g++){
            if(f == g){continue;}
            if(f != g){
                int act_sz_other = network[g]->get_nbeads();  
                for(int j = 0; j < act_sz_other; j++){
                    x1 = network[f]->get_bead(i)->get_xcm(); 
                    y1 = network[f]->get_bead(i)->get_ycm();
                    x2 = network[g]->get_bead(j)->get_xcm(); 
                    y2 = network[g]->get_bead(j)->get_ycm(); 

                    dx = x1 - x2; 
                    dy = y1 - y2; 

                    r = dist_bc(BC, dx, dy, fov[0], fov[1], delrx); 	
                    if(r == 0) { continue; } 
                    if(r <= rmax){
                        Fx1 = 2*dx*a*b*((1/r)-b); 
                        Fx2 = -Fx1; 
                        Fy1 = 2*dy*a*b*((1/r)-b); 
                        Fy2 = -Fy1; 

                        network[f]->update_forces(i,Fx1,Fy1); 
                        network[g]->update_forces(j,Fx2,Fy2);  
                    }
                    else{
                        Fx1 = 0; 
                        Fx2 = 0;
                        Fy1 = 0;
                        Fy2 = 0;  

                        network[f]->update_forces(i,Fx1,Fy1);
                        network[g]->update_forces(j,Fx2,Fy2);
                    }
                }
            } 
        }
    }
}

vector<vector<double> > filament_ensemble::spring_spring_intersections(double len, double prob){

    vector< vector<double> > itrs;
    array<double, 2> r1, r2, s1, s2, direc;
    boost::optional<array<double, 2> > inter;
    string bcf1; 
    for (unsigned int f1 = 0; f1 < network.size(); f1++){
        
        for (int l1 = 0; l1 < network[f1]->get_nsprings(); l1++){

            r1 = {{network[f1]->get_spring(l1)->get_hx()[0], network[f1]->get_spring(l1)->get_hy()[0]}};
            r2 = {{network[f1]->get_spring(l1)->get_hx()[1], network[f1]->get_spring(l1)->get_hy()[1]}};
            bcf1 = network[f1]->get_BC();
            for (unsigned int f2 = f1+1; f2 < network.size(); f2++){
                
                for (int l2 = 0; l2 < network[f2]->get_nsprings(); l2++){

                    if (f1 == f2 && fabs(double(l1) - double(l2)) < 2){ //springs should be at least two away to get crosslinked
                        continue;
                    }

                    s1 = {{network[f2]->get_spring(l2)->get_hx()[0], network[f2]->get_spring(l2)->get_hy()[0]}};
                    s2 = {{network[f2]->get_spring(l2)->get_hx()[1], network[f2]->get_spring(l2)->get_hy()[1]}};

                    inter = seg_seg_intersection_bc(bcf1, delrx, fov, r1, r2, s1, s2);
                    
                    if (inter && rng(0,1) <= prob){
                        direc = network[f2]->get_spring(l2)->get_direction();
                        itrs.push_back({inter->at(0), inter->at(1), len*direc[0], len*direc[1], 
                                double(f1), double(f2), double(l1), double(l2)}); 
                    }
                }
            }
        }
    }
    return itrs;
}
////////////////////////////////////////
///SPECIFIC FILAMENT IMPLEMENTATIONS////
////////////////////////////////////////

filament_ensemble::filament_ensemble(int npolymer, int nbeads_min, int nbeads_extra, double nbeads_extra_prob, 
        array<double,2> myfov, array<int,2> mynq, double delta_t, double temp,
        double rad, double vis, double spring_len, vector<array<double, 3> > pos_sets, double stretching, double ext, double bending, 
        double frac_force, string bc, double seed) {
    
    fov = myfov;
    view[0] = 1;//(fov[0] - 2*nbeads*len)/fov[0];
    view[1] = 1;//(fov[1] - 2*nbeads*len)/fov[1];
    nq = mynq;
    half_nq = {{nq[0]/2, nq[1]/2}};
    
    double nbeads_mean = nbeads_min + nbeads_extra*nbeads_extra_prob;
    
    visc=vis;
    spring_rest_len = spring_len;
    dt = delta_t;
    temperature = temp;
    shear_stop = 1e10;
    shear_dt = dt;
    t = 0;
    delrx = 0;
    
    if (seed == -1){
        straight_filaments = true;
    }else{
        srand(seed);
    }
    
    
    cout<<"DEBUG: Number of filament:"<<npolymer<<"\n";
    cout<<"DEBUG: Avg number of monomers per filament:"<<nbeads_mean<<"\n"; 
    cout<<"DEBUG: Monomer Length:"<<rad<<"\n"; 
   
    int nbeads = 0;
    binomial_distribution<int> distribution(nbeads_extra, nbeads_extra_prob);
    default_random_engine generator(seed+2);
    nsprings_per_fil_max = 0;

    int s = pos_sets.size();
    double x0, y0, phi0;
    for (int i=0; i<npolymer; i++) {
        if ( i < s ){
            network.push_back(new filament(pos_sets[i], nbeads, fov, nq,
                        visc, dt, temp, straight_filaments, rad, spring_rest_len, stretching, ext, bending, frac_force, bc) );
        }else{
            x0 = rng(-0.5*(view[0]*fov[0]),0.5*(view[0]*fov[0])); 
            y0 = rng(-0.5*(view[1]*fov[1]),0.5*(view[1]*fov[1]));
            phi0 =  rng(0, 2*pi);
            
            nbeads = nbeads_min + distribution(generator);
            network.push_back(new filament({{x0,y0,phi0}}, nbeads, fov, nq, visc, dt, temp, straight_filaments, rad, spring_rest_len, stretching, ext, bending, frac_force, bc) );
            nsprings_per_fil_max = max(nsprings_per_fil_max, nbeads - 1);
        }
    }
    
    //Neighbor List Initialization
    quad_off_flag = false;
    
    //this->nlist_init();
    this->nlist_init_serial();
    
    pe_stretch = 0;
    pe_bend = 0;
    
    fls = { };
}

filament_ensemble::filament_ensemble(double density, array<double,2> myfov, array<int,2> mynq, double delta_t, double temp,
        double rad, double vis, int nbeads, double spring_len, vector<array<double, 3> > pos_sets, double stretching, double ext, double bending, 
        double frac_force, string bc, double seed, double RMAX, double A) {
    
    fov = myfov;
    view[0] = 1;//(fov[0] - 2*nbeads*len)/fov[0];
    view[1] = 1;//(fov[1] - 2*nbeads*len)/fov[1];
    nq = mynq;
 
    half_nq = {{nq[0]/2, nq[1]/2}};
    BC = bc; 
    
    visc=vis;
    spring_rest_len =spring_len;
    int npolymer=int(ceil(density*fov[0]*fov[1]) / nbeads);
    dt = delta_t;
    temperature = temp;
    shear_stop = 1e10;
    shear_dt = dt;
    t = 0;
    delrx = 0;
    rmax = RMAX;
    kexv = A; 
 
    if (seed == -1){
        straight_filaments = true;
    }else{
        srand(seed);
    }
    
    
    cout<<"DEBUG: Number of filament:"<<npolymer<<"\n";
    cout<<"DEBUG: Number of monomers per filament:"<<nbeads<<"\n"; 
    cout<<"DEBUG: Monomer Length:"<<rad<<"\n"; 
   

    int s = pos_sets.size();
    double x0, y0, phi0;
    for (int i=0; i<npolymer; i++) {
        if ( i < s ){
            network.push_back(new filament(pos_sets[i], nbeads, fov, nq,
                        visc, dt, temp, straight_filaments, rad, spring_rest_len, stretching, ext, bending, frac_force, bc) );
        }else{
            x0 = rng(-0.5*(view[0]*fov[0]),0.5*(view[0]*fov[0])); 
            y0 = rng(-0.5*(view[1]*fov[1]),0.5*(view[1]*fov[1]));
            phi0 =  rng(0, 2*pi);
            //phi0=atan2(1+x0-y0*y0, -1-x0*x0+y0); // this is just the first example in mathematica's streamplot documentation
            network.push_back(new filament({{x0,y0,phi0}}, nbeads, fov, nq, visc, dt, temp, straight_filaments, rad, spring_rest_len, stretching, ext, bending, frac_force, bc) );
        }
    }
    
    //Neighbor List Initialization
    quad_off_flag = false;
    nsprings_per_fil_max = nbeads-1;

    //this->nlist_init();
    this->nlist_init_serial();
    
    pe_stretch = 0;
    pe_bend = 0;
    pe_exv = 0; 
    ke_vir = 0; 
    //ke_exv = 0;
    //N = 0;  
    
    fls = { };
}

filament_ensemble::filament_ensemble(vector<vector<double> > beads, array<double,2> myfov, array<int,2> mynq, double delta_t, double temp,
        double vis, double spring_len, double stretching, double ext, double bending, double frac_force, string bc, double RMAX, double A) {
    
    fov = myfov;

    BC = bc; 

    visc=vis;
    spring_rest_len = spring_len;
    dt = delta_t;
    temperature = temp;
    t = 0;
    delrx = 0;
    rmax = RMAX; 
    kexv = A; 
    view[0] = 1;
    view[1] = 1;

    int s = beads.size(), sa, j;
    int fil_idx = 0;
    vector<bead *> avec;
    
    nq = mynq;
    nsprings_per_fil_max = 0;
    for (int i=0; i < s; i++){
        
        if (beads[i][3] != fil_idx && avec.size() > 0){
            
            network.push_back( new filament( avec, fov, nq, spring_rest_len, stretching, ext, bending, delta_t, temp, frac_force, 0, bc) );
            
            sa = avec.size();
            nsprings_per_fil_max = max(nsprings_per_fil_max, sa - 1);
            for (j = 0; j < sa; j++) delete avec[j];
            avec.clear();
            
            fil_idx = beads[i][3];
        }
        avec.push_back(new bead(beads[i][0], beads[i][1], beads[i][2], vis));
    }

    sa = avec.size();
    if (sa > 0)
        network.push_back( new filament( avec, fov, nq, spring_rest_len, stretching, ext, bending, delta_t, temp, frac_force, 0, bc) );
    
    for (j = 0; j < sa; j++) delete avec[j];
    avec.clear();
   
    quad_off_flag = false;
    //this->nlist_init();
    this->nlist_init_serial();
    this->update_energies();

    pe_stretch = 0;
    pe_bend = 0;
    pe_exv = 0;
    ke_vir = 0; 
    //ke_exv = 0;
    //N = 0;    

    fls = { };
}

void filament_ensemble::set_growing(double kgrow, double lgrow, double l0min, double l0max, int nsprings_max)
{
    nsprings_per_fil_max = nsprings_max;
    for (int i = 0; i < int(network.size()); i++){
        network[i]->set_kgrow(kgrow);
        network[i]->set_lgrow(lgrow);
        network[i]->set_l0_min(l0min);
        network[i]->set_l0_max(l0max);
        network[i]->set_nsprings_max(nsprings_max);
    }
}
