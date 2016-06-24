/*
 *  filament_ensemble.cpp
 *  
 *
 *  Created by Simon Freedman and Shiladitya Banerjee
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
    
    for (int x = 0; x < nq[0]; x++){
        for (int y = 0; y < nq[1]; y++){
            delete links_per_quad[x]->at(y);
        }
        delete links_per_quad[x];
        delete n_links_per_quad[x];
    }
    
    for (int i = 0; i < s; i++){
        delete network[i];
    }
    
}

template <class filament_type>
vector<filament_type *>* filament_ensemble<filament_type>::get_network()
{
    return &network;
}

template <class filament_type>
filament_type * filament_ensemble<filament_type>::get_filament(int index)
{
    return network[index];
}

template<class filament_type>
void filament_ensemble<filament_type>::turn_quads_off()
{
    quad_off_flag = true;
}

template<class filament_type>
void filament_ensemble<filament_type>::nlist_init_serial()
{
    for (int x = 0; x < nq[0]; x++){
        links_per_quad.push_back(new vector< vector<array<int, 2> >* >(nq[1]+1));   
        n_links_per_quad.push_back(new vector<int>(nq[1]+1));
        for (int y = 0; y < nq[1]; y++){
            links_per_quad[x]->at(y) = new vector<array<int, 2> >(max_links_per_quad);
            n_links_per_quad[x]->at(y) = 0;
        }
    }
}

template <class filament_type> 
void filament_ensemble<filament_type>::quad_update_serial()
{
    int n_quads, net_sz = int(network.size());
    vector<vector<array<int, 2> > > q;
    int x, y;
    for (x = 0; x < nq[0]; x++)
        for (y = 0; y < nq[1]; y++)
            n_links_per_quad[x]->at(y) = 0;
    
    for (int f = 0; f < net_sz; f++){
        q = network[f]->get_quadrants();
        for (int l = 0; l < network[f]->get_nlinks(); l++){
            n_quads = int(q[l].size());
            for (int i = 0; i < n_quads; i++){
                x = q[l][i][0];
                y = q[l][i][1];
                links_per_quad[x]->at(y)->at( n_links_per_quad[x]->at(y) ) = {f,l};
                n_links_per_quad[x]->at(y)++;
            }
        }
    }

}

template <class filament_type>
void filament_ensemble<filament_type>::update_dist_map(map<array<int,2>, double>& t_map, const array<int, 2>& mq, double x, double y){
    
    array<int, 2> fl;
    double dist;
    if(n_links_per_quad[mq[0]]->at(mq[1]) != 0 ){
        
        for (unsigned int i = 0; i < n_links_per_quad[mq[0]]->at(mq[1]); i++){

            fl = links_per_quad[mq[0]]->at(mq[1])->at(i);
            dist = network[fl[0]]->get_link(fl[1])->get_distance(network[fl[0]]->get_BC(), delrx, x, y);
            
            if (t_map.find(fl) == t_map.end() || t_map[fl] > dist)
                t_map[fl] = dist;
        }
    }

}

//given motor head position, return a map between  
//  the INDICES (i.e., {i, j} where i is the filament index and j is the link index)
//  and their corresponding DISTANCES to the link at that distance 
template <class filament_type>
map<array<int,2>,double> filament_ensemble<filament_type>::get_dist(double x, double y)
{
    map<array<int, 2>, double> t_map;
    int mqx = min(coord2quad_floor(fov[0], nq[0], x), nq[0] - 1);
    int mqy = min(coord2quad_floor(fov[1], nq[1], y), nq[1] - 1);
    update_dist_map(t_map, {mqx, mqy}, x, y);
    if (mqx + 1 < nq[0])
        update_dist_map(t_map, {mqx + 1, mqy}, x, y);
    if (mqy + 1 < nq[1])
        update_dist_map(t_map, {mqx, mqy + 1}, x, y);
    if (mqx + 1 < nq[0] && mqy + 1 < nq[1])
        update_dist_map(t_map, {mqx + 1, mqy + 1}, x, y);

    return t_map;
}

template <class filament_type>
map<array<int,2>,double> filament_ensemble<filament_type>::get_dist_all(double x, double y)
{
    map<array<int, 2>, double> t_map;
    for (int f = 0; f < network.size(); f++){
        for (int l=0; l < network[f]->get_nlinks(); l++){
            t_map[{f,l}]=network[f]->get_link(l)->get_distance(network[f]->get_BC(), delrx, x, y);
        }
    }
    
    return t_map;
}

template <class filament_type>
array<double,2> filament_ensemble<filament_type>::get_intpoints(int fil, int link, double xp, double yp)
{
    return network[fil]->get_link(link)->get_intpoint(network[0]->get_BC(), delrx, xp, yp);
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
void filament_ensemble<filament_type>::update_positions()
{
    int net_sz = int(network.size());
    for (int f = 0; f < net_sz; f++)
    {
        //if (f==0) cout<<"\nDEBUG: update_positions: using "<<omp_get_num_threads()<<" cores";  
        network[f]->update_positions();
    }

}

template <class filament_type> 
void filament_ensemble<filament_type>::update_positions_range(int lo, int hi)
{
    for (unsigned int f = 0; f < network.size(); f++)
    {
        network[f]->update_positions_range(lo, hi);
    }

}

template <class filament_type> 
void filament_ensemble<filament_type>::write_actins(ofstream& fout)
{
    for (unsigned int i=0; i<network.size(); i++) {
        fout<<network[i]->write_actins(i);
    } 
}

template <class filament_type> 
void filament_ensemble<filament_type>::write_links(ofstream& fout)
{
    for (unsigned int i=0; i<network.size(); i++) {
        fout<<network[i]->write_links(i);
    } 
}

template <class filament_type> 
void filament_ensemble<filament_type>::write_thermo(ofstream& fout){
    for (unsigned int f = 0; f < network.size(); f++)
        fout<<network[f]->write_thermo(f);
    
}

template <class filament_type> 
void filament_ensemble<filament_type>::set_shear_rate(double g)
{
    if (network.size() > 0)
        if (network[0]->get_nactins() > 0)
            shear_speed = g*fov[1] / (2*network[0]->get_actin(0)->get_friction());

    for (unsigned int f = 0; f < network.size(); f++)
    {
        network[f]->set_shear(g);
    }
}

template <class filament_type> 
void filament_ensemble<filament_type>::set_y_thresh(double y)
{
    for (unsigned int f = 0; f < network.size(); f++) network[f]->set_y_thresh(y);
}

template <class filament_type> 
void filament_ensemble<filament_type>::update_delrx(double drx)
{
    //cout<<"\nDEBUG: SHEARING"; 
    delrx = drx;
    for (unsigned int f = 0; f < network.size(); f++)
    {
        network[f]->update_delrx(drx);
    }
}

template <class filament_type> 
void filament_ensemble<filament_type>::update_d_strain(double g)
{
    //cout<<"\nDEBUG: SHEARING"; 
    for (unsigned int f = 0; f < network.size(); f++)
    {
        network[f]->update_d_strain(g);
    }
}

template <class filament_type> 
void filament_ensemble<filament_type>::update_shear()
{
    //cout<<"\nDEBUG: SHEARING"; 
    for (unsigned int f = 0; f < network.size(); f++)
    {
        network[f]->update_shear(t);
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
void filament_ensemble<filament_type>::update_energies(){
    pe_stretch = 0;
    pe_bend = 0;
    ke = 0;
    for (unsigned int f = 0; f < network.size(); f++)
    {
        ke += network[f]->get_kinetic_energy();
        pe_bend += network[f]->get_bending_energy();
        pe_stretch += network[f]->get_stretching_energy();
    }
}

template <class filament_type> 
double filament_ensemble<filament_type>::get_stretching_energy(){
    return pe_stretch;
}

template <class filament_type> 
double filament_ensemble<filament_type>::get_bending_energy(){
    return pe_bend;
}

template <class filament_type> 
void filament_ensemble<filament_type>::print_network_thermo(){
    cout<<"\nAll Fs\t:\tKE = "<<ke<<"\tPEs = "<<pe_stretch<<"\tPEb = "<<pe_bend<<"\tTE = "<<(ke+pe_stretch+pe_bend);
}

template <class filament_type> 
void filament_ensemble<filament_type>::print_filament_lengths(){
    for (unsigned int f = 0; f < network.size(); f++)
    {
        cout<<"\nF"<<f<<" : "<<network[f]->get_end2end()<<" um";
    }
}


template <class filament_type> 
bool filament_ensemble<filament_type>::is_polymer_start(int fil, int actin){

    return !(actin);

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
void filament_ensemble<filament_type>::update_forces(int f_index, int a_index, double f1, double f2){
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

template <class filament_type> 
int filament_ensemble<filament_type>::get_nactins(){
    int tot = 0;
    for (unsigned int f = 0; f < network.size(); f++)
        tot += network[f]->get_nactins();
    return tot;
}

template <class filament_type> 
int filament_ensemble<filament_type>::get_nlinks(){
    return this->get_nactins() - network.size();
}

template <class filament_type> 
int filament_ensemble<filament_type>::get_nfilaments(){
    return network.size();
}

template <class filament_type> 
double filament_ensemble<filament_type>::get_delrx(){
    return delrx;
}

template <class filament_type> 
double filament_ensemble<filament_type>::get_actin_friction(){
    
    if (network.size() > 0)
        if (network[0]->get_nactins() > 0)
            return network[0]->get_actin(0)->get_friction();
    
    return 0;
}

// Update bending forces between monomers
template <class filament_type>
void filament_ensemble<filament_type>::update_bending()
{
    int net_sz = int(network.size());
    
    for (int f = 0; f < net_sz; f++)
    {
        //if (f==0) cout<<"\nDEBUG: update_bending: using "<<omp_get_num_threads()<<" cores";  
        network[f]->update_bending(t);
    }
}

template <class filament_type>
void filament_ensemble<filament_type>::update_stretching(){
    
//    vector<filament_type *> newfilaments;
    int s = network.size(); //keep it to one fracture per filament per timestep, or things get messy
    for (int f = 0; f < s; f++)
    {
        //if (f==0) cout<<"\nDEBUG: update_stretching: using "<<omp_get_num_threads()<<" cores";  
        this->update_filament_stretching(f);
    }
}

template <class filament_type>
void filament_ensemble<filament_type>::update_filament_stretching(int f){
    vector<filament_type *> newfilaments = network[f]->update_stretching(t);

    if (newfilaments.size() > 0){ //fracture event occured

        cout<<"\n\tDEBUG: fracturing filament : "<<f;
        filament * broken = network[f];     //store a pointer to the broken filament to delete it with
        network[f] = newfilaments[0];       //replace that pointer with one of the new filaments

        if (newfilaments.size() == 2) network.push_back(newfilaments[1]); //add the second filament to the top of the stack

        broken_filaments.push_back(f);      // record the index, for automatic motor detachment
        delete broken;                      // delete the old filament

    }
}

template <class filament_type>
void filament_ensemble<filament_type>::set_shear_stop(double stopT){
    shear_stop = stopT; 
}

template <class filament_type>
void filament_ensemble<filament_type>::set_shear_dt(double delT){
    shear_dt = delT; 
}

template <class filament_type>
void filament_ensemble<filament_type>::update_int_forces()
{
    this->update_stretching();
    this->update_bending();
}

/* Overdamped Langevin Dynamics Integrator (Leimkuhler, 2013) */
template <class filament_type>
void filament_ensemble<filament_type>::update()

{      
    int net_sz = network.size();
    // #pragma omp parallel for
    
    for (int f = 0; f < net_sz; f++){
      //  if (f==0) cout<<"\nDEBUG: filament updates using "<<omp_get_num_threads()<<" cores";  
        this->update_filament_stretching(f);
        network[f]->update_bending(t);
        network[f]->update_positions();
    }
    
    if (!quad_off_flag)
        this->quad_update_serial();
    
    this->update_energies();
    
    t += dt;

}

template<class filament_type>
vector<vector<double> > filament_ensemble<filament_type>::link_link_intersections(double len, double prob){

    vector< vector<double> > itrs;
    double ang;
    array<double, 2> r1, r2, s1, s2;
    pair<double, double> mmx1, mmy1, mmx2, mmy2;
    boost::optional<array<double, 2> > inter;
    string bcf1; 
    for (unsigned int f1 = 0; f1 < network.size(); f1++){
        
        for (unsigned int l1 = 0; l1 < network[f1]->get_nlinks(); l1++){

            r1 = {network[f1]->get_link(l1)->get_hx()[0], network[f1]->get_link(l1)->get_hy()[0]};
            r2 = {network[f1]->get_link(l1)->get_hx()[1], network[f1]->get_link(l1)->get_hy()[1]};
            bcf1 = network[f1]->get_BC();
            for (unsigned int f2 = f1+1; f2 < network.size(); f2++){
                
                for (unsigned int l2 = 0; l2 < network[f2]->get_nlinks(); l2++){

                    if (f1 == f2 && fabs(double(l1) - double(l2)) < 2){ //links should be at least two away to get crosslinked
                        continue;
                    }

                    s1 = {network[f2]->get_link(l2)->get_hx()[0], network[f2]->get_link(l2)->get_hy()[0]};
                    s2 = {network[f2]->get_link(l2)->get_hx()[1], network[f2]->get_link(l2)->get_hy()[1]};

                    inter = seg_seg_intersection_bc(bcf1, delrx, fov, r1, r2, s1, s2);
                    
                    if (inter && rng(0,1) <= prob){
                        ang = network[f2]->get_link(l2)->get_angle();
                        itrs.push_back({inter->at(0), inter->at(1), len*cos(ang), len*sin(ang), 
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

ATfilament_ensemble::ATfilament_ensemble(double density, array<double,2> myfov, array<int,2> mynq, double delta_t, double temp,
        double rad, double vis, int nactins, double link_len, vector<array<double, 3> > pos_sets, double stretching, double ext, double bending, 
        double frac_force, string bc, double seed) {
    
    fov = myfov;
    view[0] = 1;//(fov[0] - 2*nactins*link_len)/fov[0];
    view[1] = 1;//(fov[1] - 2*nactins*link_len)/fov[1];
    nq = mynq;
    half_nq = {nq[0]/2, nq[1]/2};
    
    visc=vis;
    link_ld = link_len;
    int npolymer=int(ceil(density*fov[0]*fov[1]) / nactins);
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
    cout<<"DEBUG: Number of monomers per filament:"<<nactins<<"\n"; 
    cout<<"DEBUG: Monomer Length:"<<rad<<"\n"; 
   

    int s = pos_sets.size();
    double x0, y0, phi0;
    for (int i=0; i<npolymer; i++) {
        if ( i < s ){
            network.push_back(new filament(pos_sets[i], nactins, fov, nq,
                        visc, dt, temp, straight_filaments, rad, link_ld, stretching, ext, bending, frac_force, bc) );
        }else{
            x0 = rng(-0.5*(view[0]*fov[0]),0.5*(view[0]*fov[0])); 
            y0 = rng(-0.5*(view[1]*fov[1]),0.5*(view[1]*fov[1]));
            phi0 =  rng(0, 2*pi);
            //phi0=atan2(1+x0-y0*y0, -1-x0*x0+y0); // this is just the first example in mathematica's streamplot documentation
            network.push_back(new filament({x0,y0,phi0}, nactins, fov, nq, visc, dt, temp, straight_filaments, rad, link_ld, stretching, ext, bending, frac_force, bc) );
        }
    }
    
    //Neighbor List Initialization
    quad_off_flag = false;
    max_links_per_quad              = npolymer*(nactins-1);
    max_links_per_quad_per_filament = nactins - 1;
    
    //this->nlist_init();
    this->nlist_init_serial();
    
    pe_stretch = 0;
    pe_bend = 0;
    ke = 0;

}

ATfilament_ensemble::ATfilament_ensemble(vector<vector<double> > actins, array<double,2> myfov, array<int,2> mynq, double delta_t, double temp,
        double vis, double link_len, double stretching, double ext, double bending, double frac_force, string bc) {
    
    fov = myfov;

    visc=vis;
    link_ld = link_len;
    dt = delta_t;
    temperature = temp;
    t = 0;
    delrx = 0;

    view[0] = 1;
    view[1] = 1;

    int s = actins.size(), sa, j;
    int fil_idx = 0;
    vector<actin *> avec;
    
    nq = mynq;
    
    for (int i=0; i < s; i++){
        
        if (actins[i][3] != fil_idx && avec.size() > 0){
            
            network.push_back( new filament( avec, fov, nq, link_len, stretching, ext, bending, delta_t, temp, frac_force, 0, bc) );
            
            sa = avec.size();
            for (j = 0; j < sa; j++) delete avec[j];
            avec.clear();
            
            fil_idx = actins[i][3];
        }
        avec.push_back(new actin(actins[i][0], actins[i][1], actins[i][2], vis));
    }

    sa = avec.size();
    if (sa > 0)
        network.push_back( new filament( avec, fov, nq, link_len, stretching, ext, bending, delta_t, temp, frac_force, 0, bc) );
    
    for (j = 0; j < sa; j++) delete avec[j];
    avec.clear();
   
    quad_off_flag = false;
    max_links_per_quad              = actins.size();
    max_links_per_quad_per_filament = int(ceil(actins.size() / (fil_idx + 1)))- 1;
    //this->nlist_init();
    this->nlist_init_serial();
} 

template class filament_ensemble<filament>;
