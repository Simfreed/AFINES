/*
 *  actin.cpp
 *  
 *
 *  Created by Shiladitya Banerjee on 9/3/13.
 *  Modified by Simon Freedman 9/2014
n*  Copyright 2013 University of Chicago. All rights reserved.
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
}

template <class filament_type> 
void filament_ensemble<filament_type>::quad_update()
{
    quad_fils_mm.clear();
    int qs, net_sz = int(network.size());
    vector<array<int, 2> > q;
    vector<multimap<array<int, 2>, array<int, 2> >* > fil_quads(net_sz);
//    <multimap<array<int, 2>, array<int, 2> > quad_fils_mmi;
    
    #pragma omp parallel for
    for (int f = 0; f < net_sz; f++){
        multimap<array<int, 2>, array<int, 2> > *mm = new multimap<array<int, 2>, array<int, 2> >;
        for (int l = 0; l < network[f]->get_nlinks(); l++){
            network[f]->get_link(l)->quad_update(network[f]->get_BC(), delrx); 
            q = network[f]->get_link(l)->get_quadrants();
            qs = int(q.size());
            for (int i = 0; i < qs; i++){
//                quad_fils_mm.insert(pair<array<int, 2>, array<int, 2> >(q[i], {f,l}));
                mm->insert(pair<array<int, 2>, array<int, 2> >(q[i], {f,l}));
            }
        }
        fil_quads[f] = mm;
    }

    for (int f = 0; f < net_sz; f++){
        quad_fils_mm.insert(fil_quads[f]->begin(), fil_quads[f]->end());
        delete fil_quads[f];
    }

}

/*
template <class filament_type> 
void filament_ensemble<filament_type>::quad_update()
{
    quad_fils.clear();
    vector< vector< array<int, 2> > > filament_quads;
    const int net_sz = int(network.size());
    map<int, vector< vector< array<int, 2> > > > all_filament_quads;
//    #pragma omp parallel for
    
    for (int i=0; i < net_sz; i++) { //Loop over filaments
        
        // filament_quads = network[i]->get_quadrants(); 
        all_filament_quads[i] = network[i]->get_quadrants();
        
    }
    
    for (map<int, vector< vector< array<int, 2> > > >::iterator it = all_filament_quads.begin(); it != all_filament_quads.end(); ++it){
        filament_quads = it->second;
        for (unsigned int j=0; j < filament_quads.size(); j++){ //Loop over links

            for (unsigned int k = 0; k  < filament_quads[j].size(); k++){ //Loop over quadrants of a Link   

                quad_fils[ filament_quads[j][k] ].push_back({it->first, (int)j});
                //cout<<"\nDEBUG: quad_fils[{"<<filament_quads[j][k][0]<<" , "<<filament_quads[j][k][1]<<"}] = {"<<i<<" , "<<j<<"}";
            }
        }   
    }

}*/

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

//given motor head position, return a map between  
//  the INDICES (i.e., {i, j} where i is the filament index and j is the link index)
//  and their corresponding DISTANCES to the link at that distance 
template <class filament_type>
map<array<int,2>,double> filament_ensemble<filament_type>::get_dist(double x, double y)
{
    map<array<int, 2>, double> t_map;
    int mqx = int(floor(x/fov[0]*nq[0]));
    int mqy = int(floor(y/fov[1]*nq[1]));
    update_dist_map(t_map, {mqx, mqy}, x, y);
    update_dist_map(t_map, {mqx, mqy + 1}, x, y);
    update_dist_map(t_map, {mqx + 1, mqy}, x, y);
    update_dist_map(t_map, {mqx + 1, mqy + 1}, x, y);
    
    return t_map;
}
/*
template <class filament_type>
void filament_ensemble<filament_type>::update_dist_map(map<array<int,2>, double>& t_map, const array<int, 2>& mquad, double x, double y){
    
    array<int, 2> lnk_idx;
    double dist;
    
    if(quad_fils.find(mquad) != quad_fils.end()){
        
        for (unsigned int j = 0; j < quad_fils[mquad].size(); j++){
        
            lnk_idx = quad_fils[mquad][j];
            dist = network[lnk_idx[0]]->get_link(lnk_idx[1])->get_distance(network[lnk_idx[0]]->get_BC(), delrx, x, y);
            
            if (t_map.count(lnk_idx) == 0 || t_map[lnk_idx] > dist)
                
                t_map[lnk_idx] = dist;
        }
    }

}
*/
template <class filament_type>
void filament_ensemble<filament_type>::update_dist_map(map<array<int,2>, double>& t_map, const array<int, 2>& mquad, double x, double y){
    
    array<int, 2> fl_idx;
    double dist;
    multimap<array<int, 2>, array<int, 2> >::iterator it;
    pair<multimap<array<int, 2>, array<int, 2> >::iterator, multimap<array<int, 2>, array<int, 2> >::iterator > ii;
    ii = quad_fils_mm.equal_range(mquad);
    
    for(it = ii.first; it != ii.second; ++it){
        fl_idx = it->second;
        dist = network[fl_idx[0]]->get_link(fl_idx[1])->get_distance(network[fl_idx[0]]->get_BC(), delrx, x, y);
        
        if (t_map.find(fl_idx) == t_map.end() || t_map[fl_idx] > dist) t_map[fl_idx] = dist;
    }

}
    
template <class filament_type>
array<double,2> filament_ensemble<filament_type>::get_intpoints(int fil, int link, double xp, double yp)
{
    return network[fil]->get_link(link)->get_intpoint(network[0]->get_BC(), delrx, xp, yp);
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
void filament_ensemble<filament_type>::update_positions()
{
    int net_sz = int(network.size());
    #pragma omp parallel for 
    for (int f = 0; f < net_sz; f++)
    {
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
void filament_ensemble<filament_type>::print_network_thermo(){
    double KE=0, PEs=0, PEb=0, TE=0;
    for (unsigned int f = 0; f < network.size(); f++)
    {
        KE += network[f]->get_kinetic_energy();
        //PE += network[f]->get_potential_energy();
        PEb += network[f]->get_bending_energy();
        PEs += network[f]->get_stretching_energy();
        TE += network[f]->get_total_energy();
    }
    cout<<"\nAll Fs\t:\tKE = "<<KE<<"\tPEs = "<<PEs<<"\tPEb = "<<PEb<<"\tTE = "<<TE;
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
    #pragma omp parallel for
    
    for (int f = 0; f < net_sz; f++)
    {
        network[f]->update_bending(t);
    }
}

template <class filament_type>
void filament_ensemble<filament_type>::update_stretching(){
    
//    vector<filament_type *> newfilaments;
    int s = network.size(); //keep it to one fracture per filament per timestep, or things get messy
    #pragma omp parallel for
    for (int f = 0; f < s; f++)
    {
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
void filament_ensemble<filament_type>::update(){
    
    this->update_int_forces();
    this->update_positions();
    this->quad_update();
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
    nq = mynq;

    view[0] = 1;//(fov[0] - 2*nactins*link_len)/fov[0];
    view[1] = 1;//(fov[1] - 2*nactins*link_len)/fov[1];

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
    
    max_links_per_gp =  10*nq[0]*nq[1]/(npolymer*nactin);

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
}

ATfilament_ensemble::ATfilament_ensemble(vector<vector<double> > actins, array<double,2> myfov, array<int,2> mynq, double delta_t, double temp,
        double vis, double link_len, double stretching, double ext, double bending, double frac_force, string bc) {
    
    fov = myfov;
    nq = mynq;
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
} 

template class filament_ensemble<filament>;
