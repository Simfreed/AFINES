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
    int qs = all_quads.size();    
    for (unsigned int i = 0; i < qs; i++){
        /*for (unsigned int f = 0; f < network.size(); f++){
            delete links_per_quad_per_filament[f]->at(*(all_quads[i]));
        }*/
        delete links_per_quad[*(all_quads[i])];
        delete all_quads[i];
    }
    
 /*   for (unsigned int f = 0; f < s; f++){
        delete links_per_quad_per_filament[f];
        delete n_links_per_quad_per_filament[f];
    }
   */ 
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
void filament_ensemble<filament_type>::nlist_init()
{
    array<int, 2> xy;
    for (int x = -nq[0]/2; x <= nq[0]/2; x++)
        for (int y = -nq[1]/2; y <= nq[1]/2; y++){
            all_quads.push_back(new array<int,2>{{x,y}});
            xy = *(all_quads[all_quads.size()-1]);
        }

    for (array<int, 2>* xyptr : all_quads){ 
        links_per_quad[*xyptr] = new vector<array<int,2> >(max_links_per_quad);
        n_links_per_quad[*xyptr] = 0;
    }
    
    for (unsigned int f = 0; f < network.size(); f++){
        //links_per_quad_per_filament.push_back(new unordered_map<array<int,2>, vector<int>*, boost::hash<array<int, 2> > >);
        //n_links_per_quad_per_filament.push_back(new unordered_map<array<int, 2>, int, boost::hash<array<int, 2> > >);
        links_per_quad_per_filament.push_back(new map<array<int,2>, vector<int>*>);
        n_links_per_quad_per_filament.push_back(new map<array<int, 2>, int>);
        for (array<int, 2>* xyptr : all_quads){ 
            links_per_quad_per_filament[f]->insert(
                    pair<array<int,2>, vector<int>* >(*xyptr, new vector<int>(max_links_per_quad_per_filament) ) );
            n_links_per_quad_per_filament[f]->insert(pair<array<int, 2>, int>(*xyptr, 0) );
//            cout<<"\nDEBUG: {nlinks_per_quad_per_filament["<<f<<"]->at{"<<xy[0]<<","<<xy[1]<<"} => "<<n_links_per_quad_per_filament[f]->at(xy);
        }
    }

}

template<class filament_type>
void filament_ensemble<filament_type>::nlist_init_serial()
{
    array<int, 2> xy;
    for (int x = -nq[0]/2; x <= nq[0]/2; x++)
        for (int y = -nq[1]/2; y <= nq[1]/2; y++){
            all_quads.push_back(new array<int,2>{{x,y}});
            xy = *(all_quads[all_quads.size()-1]);
        }

    for (array<int, 2>* xyptr : all_quads){ 
        links_per_quad[*xyptr] = new vector<array<int,2> >(max_links_per_quad);
    }

}
/*
template <class filament_type> 
void filament_ensemble<filament_type>::reset_n_links(int f)
{
    for (array<int, 2>* xyptr : all_quads){ 
        n_links_per_quad_per_filament[f]->at(*xyptr)=0;
    }
}
*/
template <class filament_type> 
void filament_ensemble<filament_type>::update_quads_per_filament(int f)
{
    vector<vector<array<int, 2> > > filament_quads = network[f]->get_quadrants(); 
    for (unsigned int l = 0; l < filament_quads.size(); l++){
        for (unsigned int gp = 0; gp < filament_quads[l].size(); gp++){
            //links_per_quad[at filament f][at grid point {x,y}][at the highest untaken index]          = link_index
            links_per_quad_per_filament[f]->at(filament_quads[l][gp])->at(
                    n_links_per_quad_per_filament[f]->at(filament_quads[l][gp]) ) = l;
            n_links_per_quad_per_filament[f]->at(filament_quads[l][gp]) += 1;
            if (n_links_per_quad_per_filament[f]->at(filament_quads[l][gp]) == max_links_per_quad_per_filament)
                break;
        }
    }
}

template <class filament_type> 
void filament_ensemble<filament_type>::consolidate_quads()
{
    array<int, 2> xy;
    int nl;
    const int net_sz = int(network.size());
    for (array<int, 2>* xyptr : all_quads){
        n_links_per_quad[*xyptr] = 0;
    }
    for (int f=0; f < net_sz; f++) { //Loop over filaments
        // {x,y} --> vec ints
        for (map<array<int,2>, vector<int>*>::iterator it = links_per_quad_per_filament[f]->begin(); 
                it != links_per_quad_per_filament[f]->end(); ++it)
        {
            xy = it->first; 
            nl = n_links_per_quad_per_filament[f]->at(xy);
            if (nl == 0) continue;
            for (int i = 0; i < nl; i++){
                links_per_quad[xy]->at(i + n_links_per_quad[xy] ) = {f, it->second->at(i)};
            }
            n_links_per_quad[xy] += nl;
            n_links_per_quad_per_filament[f]->at(xy) = 0;
        }
    }

   /* 
    for (array<int, 2>* xyptr : all_quads){
        xy = *xyptr;
        n_links_per_quad[xy] = 0;
        for (int f=0; f < net_sz; f++) { //Loop over filaments
            for (int i = 0; i < n_links_per_quad_per_filament[f]->at(xy); i++){
                links_per_quad[xy]->at(i + n_links_per_quad[xy] ) = {f, links_per_quad_per_filament[f]->at(xy)->at(i)};
            }
            n_links_per_quad[xy] += n_links_per_quad_per_filament[f]->at(xy);
            n_links_per_quad_per_filament[f]->at(xy) = 0;
        }
    }
    */
}
    
template <class filament_type> 
void filament_ensemble<filament_type>::quad_update()
{
    const int net_sz = int(network.size());
    
    #pragma omp parallel for
    
    for (int f=0; f < net_sz; f++) { //Loop over filaments
        if (f==0) cout<<"\nDEBUG: quad_update: using "<<omp_get_num_threads()<<" cores";  
        update_quads_per_filament(f);
    }
    
    consolidate_quads();
}

template <class filament_type>
void filament_ensemble<filament_type>::update_dist_map(map<array<int,2>, double>& t_map, const array<int, 2>& mquad, double x, double y){
    
    array<int, 2> fl;
    double dist;
    
    if(links_per_quad.count(mquad)){
        
        for (unsigned int i = 0; i < n_links_per_quad[mquad]; i++){

            fl = links_per_quad[mquad]->at(i);
            dist = network[fl[0]]->get_link(fl[1])->get_distance(network[fl[0]]->get_BC(), delrx, x, y);
            
            if (t_map.find(fl) == t_map.end() || t_map[fl] > dist)
                t_map[fl] = dist;
        }
    }

}

template <class filament_type> 
void filament_ensemble<filament_type>::quad_update_serial()
{
    int qs, net_sz = int(network.size());
    vector<array<int, 2> > q;
    
    for (array<int, 2>* xyptr : all_quads){ 
        n_links_per_quad[*xyptr]=0;
    }
    
    for (int f = 0; f < net_sz; f++){
        for (int l = 0; l < network[f]->get_nlinks(); l++){
            q = network[f]->get_link(l)->get_quadrants();
            qs = int(q.size());
            for (int i = 0; i < qs; i++){
                links_per_quad[q[i]]->at(n_links_per_quad[q[i]]) = {f,l};
                n_links_per_quad[q[i]]++;
            }
        }
    }

}

/*
template <class filament_type>
void filament_ensemble<filament_type>::update_dist_map_serial(map<array<int,2>, double>& t_map, const array<int, 2>& mquad, double x, double y){
    
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
 */   
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
        if (f==0) cout<<"\nDEBUG: update_positions: using "<<omp_get_num_threads()<<" cores";  
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
        if (f==0) cout<<"\nDEBUG: update_bending: using "<<omp_get_num_threads()<<" cores";  
     
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
        if (f==0) cout<<"\nDEBUG: update_stretching: using "<<omp_get_num_threads()<<" cores";  
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
    //this->quad_update();
    this->quad_update_serial();
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
    max_links_per_quad              = npolymer*(nactins-1);
    max_links_per_quad_per_filament = nactins - 1;
    //this->nlist_init();
    this->nlist_init_serial();

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
   
    max_links_per_quad              = actins.size();
    max_links_per_quad_per_filament = int(ceil(actins.size() / (fil_idx + 1)))- 1;
    //this->nlist_init();
    this->nlist_init_serial();
} 

template class filament_ensemble<filament>;
