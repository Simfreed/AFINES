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
#include "actin_ensemble.h"
//actin network class

actin_ensemble::actin_ensemble(){}

actin_ensemble::actin_ensemble(double density, double fovx, double fovy, int nx, int ny, double delta_t, double temp,
        double len, double vis, int nmonomer, double link_len, std::vector<double *> pos_sets, double seed) {
    
    view[0] = (fovx - 2*nmonomer*len)/fovx;
    view[1] = (fovy - 2*nmonomer*len)/fovy;
    fov[0]=fovx;
    fov[1]=fovy;
    nq[0]=nx;
    nq[1]=ny;
    rho=density;
    visc=vis;
    ld=len;//rng_n(len,1.0);
    link_ld = link_len;
    npolymer=int(ceil(density*fov[0]*fov[1]) / nmonomer);
    dt = delta_t;
    temperature = temp;
   
    if (seed == -1){
        straight_filaments = true;
    }else{
        srand(seed);
    }

    std::cout<<"DEBUG: Number of filament:"<<npolymer<<"\n";
    std::cout<<"DEBUG: Number of monomers per filament:"<<nmonomer<<"\n"; 
    std::cout<<"DEBUG: Monomer Length:"<<ld<<"\n"; 
    
    for (int i=0; i<npolymer; i++) {
        int s = pos_sets.size();
        if ( i < s){
            this->add_polymer(pos_sets[i][0], pos_sets[i][1], pos_sets[i][2], i, nmonomer);
        }else{
            this->add_polymer(rng(-0.5*(view[0]*fov[0]-ld),0.5*(view[0]*fov[0]-ld)), rng(-0.5*(view[1]*fov[1]-ld),0.5*(view[1]*fov[1]-ld)), rng(0, 2*pi), i, nmonomer);
        }
    }
}

actin_ensemble::~actin_ensemble(){ 
    std::cout<<"DELETING ACTIN_ENSEMBLE\n";
    
    int s = network.size();
    for (int i = 0; i < s; i++){
        delete network[i];
    }
    
    network.clear();
    actin_link_map.clear();
};


void actin_ensemble::add_polymer(double startx, double starty, double phi0, int pol_index, int nmon)
{
    //the start of the polymer: 
    actin * a = new actin(startx, starty, phi0, ld, fov[0], fov[1], nq[0], nq[1], visc);
    add_monomer(a, pol_index);
    
    double  xcm, ycm, lphi, phi;
    phi = phi0;
    for (int j = 1; j < nmon; j++) {

        // Calculate the Next rod on the actin polymer--  continues from the link
        if (!straight_filaments){ //constrain phi to be <= 90 degree difference in either direction
            phi += rng(-pi/2,pi/2);
        }
            
        lphi = (phi + network.back()->get_angle())/2;
        xcm = network.back()->get_end()[0] + link_ld*cos(lphi) + ld*0.5*cos(phi);
        ycm = network.back()->get_end()[1] + link_ld*sin(lphi) + ld*0.5*sin(phi);

        // Check that this monomer is in the field of view, otherwise start a new polymer:
        if ( xcm > (0.5*(fov[0] - ld)) || xcm < (-0.5*(fov[0] - ld)) 
                || ycm > (0.5*(fov[1] - ld)) || ycm < (-0.5*(fov[1] - ld)) )
        {
            std::cout<<"DEBUG:"<<j+1<<"th monomer of "<<pol_index<<"th polymer outside field of view; stopped building polymer\n";
            break;
        }else{
            // Add the actin monomer
            a = new actin(xcm, ycm, phi, ld, fov[0], fov[1], nq[0], nq[1], visc) ;
            add_monomer(a, pol_index); 
        } 

    }
}

void actin_ensemble::add_monomer(actin * a, int polymer){
    network.push_back(a);
    mono_map[polymer].push_back(network.size() - 1);
    quad_update_monomer(network.size() - 1);
}

void actin_ensemble::quad_update_monomer(int i){
    
    std::vector<std::vector<int> > tmp_quads=network[i]->get_quadrants();
    for (unsigned int xindex=0; xindex<tmp_quads[0].size(); xindex++) {
        for (unsigned int yindex=0; yindex<tmp_quads[1].size(); yindex++) {
            quad_fils[tmp_quads[0][xindex]][tmp_quads[1][yindex]].push_back(i);
        }
    }

}
void actin_ensemble::quad_update()
{
    quad_fils.clear();
    for (unsigned int i=0; i<network.size(); i++) {
        quad_update_monomer(i);
    }
}

std::vector<actin *>* actin_ensemble::get_network()
{
    return &network;
}

//given motor quadrant, return the indices and the distances to all actin filaments in the neighboring quadrants 
std::map<int,double> actin_ensemble::get_dist(double x, double y)
{
    int nxx=int(floor(x/fov[0]*nq[0]));
    int nyy=int(floor(y/fov[1]*nq[1]));
    t_map.clear();
    if(!quad_fils[nxx][nyy].empty())
    {
        for(std::vector<int>::iterator it=quad_fils[nxx][nyy].begin(); it<quad_fils[nxx][nyy].end(); it++)
        {   
            t_map[*it] = network[*it]->get_distance(x,y);
        }
    }
    return t_map;
}

double* actin_ensemble::get_direction(int index)
{
    return network[index]->get_direction();
}


double* actin_ensemble::get_intpoints(int index, double xp, double yp)
{
    return network[index]->get_intpoint(xp,yp);
}

double actin_ensemble::get_int_direction(int index, double xp, double yp)
{
    return network[index]->get_int_angle(xp,yp);
}

double actin_ensemble::get_xcm(int index)
{
    return network[index]->get_xcm();
}

double actin_ensemble::get_ycm(int index)
{
    return network[index]->get_ycm();
}

double actin_ensemble::get_angle(int index)
{
    return network[index]->get_angle();
}

double actin_ensemble::get_alength(int index)
{
    return network[index]->get_length();
}

double * actin_ensemble::get_start(int index)
{
    return network[index]->get_start();
}

double * actin_ensemble::get_end(int index)
{
    return network[index]->get_end();
}

double * actin_ensemble::get_forces(int index)
{
    return network[index]->get_forces();
}

void actin_ensemble::set_straight_filaments(bool is_straight)
{
    straight_filaments = is_straight;
}

void actin_ensemble::update()
{
    ///Maybe change 6 to 4 for 2d
    double vpar, vperp, vx, vy, omega, alength, xnew, ynew, phinew, a_ends[4]; 
    
    for (unsigned int i=0; i<network.size(); i++) {
        
        double * fric = network[i]->get_friction();
        vpar=(network[i]->get_forces()[0])/fric[0]  + sqrt(4*temperature/(dt*fric[0]))*rng_n(0,1);
        vperp=(network[i]->get_forces()[1])/fric[1] + sqrt(4*temperature/(dt*fric[1]))*rng_n(0,1);
        vx=vpar*cos(network[i]->get_angle())-vperp*sin(network[i]->get_angle());
        vy=vpar*sin(network[i]->get_angle())+vperp*cos(network[i]->get_angle());
        omega=network[i]->get_forces()[2]/fric[2] + sqrt(4*temperature/(dt*fric[2]))*rng_n(0,1);
        delete[] fric;

        alength=network[i]->get_length();

        xnew=network[i]->get_xcm()+dt*vx;
        ynew=network[i]->get_ycm()+dt*vy;
        phinew=network[i]->get_angle()+dt*omega;


        a_ends[0]=xnew-alength*0.5*cos(phinew);
        a_ends[1]=ynew-alength*0.5*sin(phinew);
        a_ends[2]=xnew+alength*0.5*cos(phinew);
        a_ends[3]=ynew+alength*0.5*sin(phinew);

        if (a_ends[0]<=-fov[0]*0.5 || a_ends[0]>=fov[0]*0.5 || a_ends[2]<=-fov[0]*0.5 || a_ends[2]>=fov[0]*0.5)
        {
            vx=-vx;//xnew=network[i]->get_xcm()-dt*vx;//  
            omega=-omega;//phinew=network[i]->get_angle()-dt*omega;//omega=-omega;

        }
        if (a_ends[1]<=-fov[1]*0.5 || a_ends[1]>=fov[1]*0.5 || a_ends[3]<=-fov[1]*0.5 || a_ends[3]>=fov[1]*0.5)
        {
            vy=-vy;//ynew=network[i]->get_ycm()-dt*vy;
            omega=-omega;//phinew=network[i]->get_angle()-dt*omega;
        }

        xnew=network[i]->get_xcm()+dt*vx;
        ynew=network[i]->get_ycm()+dt*vy;
        phinew=network[i]->get_angle()+dt*omega;
        network[i]->set_xcm(xnew);
        network[i]->set_ycm(ynew);
        network[i]->set_phi(phinew);
        network[i]->update();
    }


}

void actin_ensemble::update_forces(int index, double f1, double f2, double f3)
{
    network[index]->update_force(f1,f2,f3);
}


void actin_ensemble::write(std::ofstream& fout)
{
    for (unsigned int i=0; i<network.size(); i++) {
        fout<<network.at(i)->write();
    } 
}
/*
   wikipedia says < cos(theta) > = Exp[-L/P]
   where L is the length of the polymer and P is the persistence length
   and theta is the angle between subsequent rods
   This function will return a map of x --> <cos theta>
   where x is a distance and <cos theta> is the correlation between cosines of the angles between subsequent
   rods. 
   */ 
std::vector<double> actin_ensemble::get_angle_correlation(int polymer_index)
{
    double sum = 0, phi0, phi1;//, phi2;
    std::vector<double> corr;
    std::vector<int> * monos = &mono_map[polymer_index]; 
    phi0 = network.at(monos->at(0))->get_angle();
    for(unsigned int i = 0; i < monos->size(); i++){

        phi1 = network.at(monos->at(i))->get_angle(); 
        //phi2 = network.at(monos->at(i+1))->get_angle(); 

        sum += cos(phi1 - phi0);

        corr.push_back(sum/(i+1));

    }
    return corr;

}

std::map<int, std::vector<double> > actin_ensemble::get_all_angle_correlations()
{
    std::map<int, std::vector<double> > corrs;

    for(unsigned int i = 0; i < mono_map.size(); i++){
        corrs[i] = this->get_angle_correlation(i);
    }

    return corrs;
}

/* Based on derivation in :
 * Flexural Rigidity of Microtubules and Actin Filaments etc. 
 * Gittes, Micky, Nettelton Howard
 * Journal of Cell Biology; Feb 1993
 */
double actin_ensemble::get_fourier_mode(int n, int polymer_index){

    std::vector<int> mons = mono_map[polymer_index];
    int s = mons.size();
    double L = s * ld, sum = 0, phi = 0, sk = 0;

    for (int i = 0; i < s; i++){

        phi = network[ mons[ i ] ]->get_angle();
        sk = i * ld + ld/2;
        sum += phi * ld * cos( n * pi * sk / L);

    }

    return sqrt(2.0/L)*sum;

}

// Adds Links between monomers
void actin_ensemble::connect_polymers(link_ensemble * links, double link_length, double
        stretching_stiffness, double bending_stiffness, std::string link_color)
{
    int mono1, mono2;
    Link * l;
    for (unsigned int i = 0; i < mono_map.size(); i++){
        
        mono1 = mono_map[i][0];
        l = new Link( link_length, stretching_stiffness, bending_stiffness, this, -1, mono1,  link_color); 
        links->add_link( l );
        actin_link_map[-1][mono1] = l;
        for (unsigned int j = 1; j < mono_map[i].size(); j++){
            mono2 = mono_map[i][j];
            l = new Link( link_length, stretching_stiffness, bending_stiffness, this, mono1, mono2, link_color); 
            links->add_link( l );
            actin_link_map[mono1][mono2] = l;
        //    std::cout<<"DEBUG: created "<<j<<"th link between monomers "<<mono1<<" and "<<mono2<<"\n";
        //    std::cout<<"DEBUG: "<<l->to_string();
            mono1 = mono2;
        }
        l = new Link( link_length, stretching_stiffness, bending_stiffness, this, mono1, -1, link_color); 
        links->add_link( l );
        actin_link_map[mono1][-1] = l;
    }

}

void actin_ensemble::update_polymer_bending(int polymer_index)
{

    Link * lft_lnk, * ctr_lnk, * rt_lnk;
    double forcex, forcey, force_par, force_perp, lft_trq, rt_trq;
    int mono1, mono2, mono3;
    std::vector<double> node_forces_x, node_forces_y;
    std::vector<int> * monomers = &mono_map[polymer_index];
    
    if (monomers->size() < 2){ //no bending energy!
        return;
    }
    
    //initialize all NODE forces to be 0
    for (unsigned int j = 0; j <= monomers->size(); j++){
        node_forces_x.push_back(0);
        node_forces_y.push_back(0);
    }
    
    //Calculate the force at each Link position as outlined by Nedelec, Foethke (2007)
    //Keep the forces at the ends of each filament 0
    
    mono1 = monomers->at(0);
    mono2 = monomers->at(1);
    lft_lnk = actin_link_map[-1][mono1];
    ctr_lnk = actin_link_map[mono1][mono2];
    
    for (unsigned int j = 2; j <= monomers->size(); j++){

        if( j == monomers->size() )
            rt_lnk = actin_link_map[mono2][-1];
        else{ 
            mono3 = monomers->at(j);
            rt_lnk = actin_link_map[mono2][mono3];
        }
        forcex =    lft_lnk->get_kb() * lft_lnk->get_posx() 
              - 2 * ctr_lnk->get_kb() * ctr_lnk->get_posx() 
              +     rt_lnk->get_kb()  * rt_lnk->get_posx(); 
        
        forcey =    lft_lnk->get_kb() * lft_lnk->get_posy() 
              - 2 * ctr_lnk->get_kb() * ctr_lnk->get_posy() 
              +     rt_lnk->get_kb()  * rt_lnk->get_posy();

        node_forces_x[j-2] += -1 * forcex;
        node_forces_x[j-1] +=  2 * forcex;
        node_forces_x[j  ] += -1 * forcex;

        node_forces_y[j-2] += -1 * forcey;
        node_forces_y[j-1] +=  2 * forcey;
        node_forces_y[j  ] += -1 * forcey;

        lft_lnk = ctr_lnk;
        ctr_lnk = rt_lnk;
        mono2 = mono3;
    }
    //Calculate the forces at each center of mass of the monomers
    //Update the monomer

    lft_lnk = actin_link_map[-1][monomers->at(0)];

    for (unsigned int j = 0; j < monomers->size(); j++){
    
        forcex = (node_forces_x[j] + node_forces_x[j+1]) / 2;
        forcey = (node_forces_y[j] + node_forces_y[j+1]) / 2;

        force_par   =  forcex*this->get_direction(monomers->at(j))[0] + forcey*this->get_direction(monomers->at(j))[1];
        force_perp  = -forcex*this->get_direction(monomers->at(j))[1] + forcey*this->get_direction(monomers->at(j))[0];
        
        if (j == 0)
            lft_trq = 0;
        else
            lft_trq = cross(this->get_xcm(monomers->at(j)) - lft_lnk->get_posx(),
                            this->get_ycm(monomers->at(j)) - lft_lnk->get_posy(), forcex, forcey);

        if (j == monomers->size()-1)
            rt_trq = 0;
        else{
            rt_lnk = actin_link_map[monomers->at(j)][monomers->at(j+1)];
            rt_trq = cross(this->get_xcm(monomers->at(j)) - rt_lnk->get_posx(),
                           this->get_ycm(monomers->at(j)) - rt_lnk->get_posy(), forcex, forcey);
        }
        
        this->update_forces(monomers->at(j), force_par, force_perp, lft_trq + rt_trq);

        lft_lnk = rt_lnk;
    }
}

// Update bending forces between monomers
void actin_ensemble::update_bending(){
    
    for (unsigned int p = 0; p < mono_map.size(); p++)
    {
        this->update_polymer_bending(p);
    }
}

void actin_ensemble::clear_actin_link_map(){
    
    actin_link_map.clear();

}
