/*
 *  actin.cpp
 *  
 *
 *  Created by Shiladitya Banerjee on 9/3/13.
 *  Copyright 2013 University of Chicago. All rights reserved.
 *
 */

#include "globals.h"
#include "actin_ensemble.h"
#include "link_ensemble.h"
#include "Link.h"
//actin network class

actin_ensemble::actin_ensemble(double density, double fovx, double fovy, int nx, int ny, double len, double vis, int nmonomer, double link_len)
{
    view=0.9;
    fov[0]=fovx;
    fov[1]=fovy;
    nq[0]=nx;
    nq[1]=ny;
    rho=density;
    av_vel=0;
    visc=vis;
    //        nmonomer_min = 100; //hard coded number of min/max monomers per filament
    //        nmonomer_max = nmonomer_min;
    //        int nmonomer = (nmonomer_max + nmonomer_min)/2;
    npolymer=int(ceil(density*fov[0]*fov[1]) / nmonomer);
    ld=len;//rng_n(len,1.0);
    link_ld = link_len;
    std::cout<<"DEBUG: Number of filament:"<<npolymer<<"\n";
    std::cout<<"DEBUG: Number of monomers per filament:"<<nmonomer<<"\n"; 
    std::cout<<"DEBUG: Monomer Length:"<<ld<<"\n"; 
    
    double  xcm, ycm, theta;
    for (int i=0; i<npolymer; i++) {
        theta=rng(0,2*pi);
        //nmonomer = (int) rng(nmonomer_min, nmonomer_max);
        //the start of the polymer: 
                    network.push_back(actin(rng(-0.5*(view*fovx-ld),0.5*(view*fovx-ld)), rng(-0.5*(view*fovy-ld),0.5*(view*fovy-ld)),
                                theta,ld,fov[0],fov[1],nq[0],nq[1],visc));
        //std::cout<<"WARNING: STARTING ACTIN FILAMENT POSITION CHOSEN DETERMINISTICALLY\n";
        //network.push_back(actin(0,0,theta,ld,fov[0],fov[1],nq[0],nq[1],visc));
        //Add the quadrants of the first rod
        std::vector<std::vector<int> > tmp_quads=network.back().get_quadrants();
        for (unsigned int xindex=0; xindex<tmp_quads[0].size(); xindex++) {
            for (unsigned int yindex=0; yindex<tmp_quads[1].size(); yindex++) {
                quad_fils[tmp_quads[0][xindex]][tmp_quads[1][yindex]].push_back(network.size()-1);
            }
        }
        // add monomers to the polymer
        mono_map[i] = empty_vector; 
        for (int j=0; j<nmonomer-1; j++) {

            // Calculate a link (represented as a dead motor) at the end of the rod
            ltheta = theta; //rng(0, 2*pi);

            // Calculate the Next rod on the actin polymer--  continues from the link
            //theta = rng(0,2*pi);
            xcm = network.back().getendpts()[2] + link_ld*cos(ltheta) + ld*0.5*cos(theta);
            ycm = network.back().getendpts()[3] + link_ld*sin(ltheta) + ld*0.5*sin(theta);

            // Check that this monomer is in the field of view, otherwise start a new polymer:
            if ( xcm > (0.5*(view*fovx - ld)) || xcm < (-0.5*(view*fovx - ld)) 
                    || ycm > (0.5*(view*fovy - ld)) || ycm < (-0.5*(view*fovy - ld)) )
            {
                std::cout<<"\nDEBUG:"<<j+1<<"th monomer of "<<i<<"th polymer outside field of view; stopped building polymer\n";
                break;
            }else{

                mono_map[ i ].push_back(network.size() -1);
                
                // Add the actin monomer
                network.push_back( actin(xcm, ycm, theta, ld, fov[0], fov[1], nq[0], nq[1], visc) );

                // Add it's quadrants:
                std::vector<std::vector<int> > tmp_quads=network.back().get_quadrants();
                for (unsigned int xindex=0; xindex<tmp_quads[0].size(); xindex++) {
                    for (unsigned int yindex=0; yindex<tmp_quads[1].size(); yindex++) {
                        quad_fils[tmp_quads[0][xindex]][tmp_quads[1][yindex]].push_back(network.size()-1);
                    }
                }
            } 

        }
    }


}

actin_ensemble::~actin_ensemble(){ };

void actin_ensemble::quad_update()
{
    quad_fils.clear();
    for (unsigned int i=0; i<network.size(); i++) {
        std::vector<std::vector<int> > tmp_quads=network[i].get_quadrants();
        for (unsigned int xindex=0; xindex<tmp_quads[0].size(); xindex++) {
            for (unsigned int yindex=0; yindex<tmp_quads[1].size(); yindex++) {
                quad_fils[tmp_quads[0][xindex]][tmp_quads[1][yindex]].push_back(i);
            }
        }
    }


}

std::vector<actin>* actin_ensemble::get_network()
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
            t_map[*it] = network[*it].get_distance(x,y);
        }
    }
    return t_map;
}

double* actin_ensemble::get_direction(int index)
{
    return network[index].get_direction();
}


double* actin_ensemble::get_intpoints(int index, double xp, double yp)
{
    return network[index].get_intpoint(xp,yp);
}

double actin_ensemble::get_int_direction(int index, double xp, double yp)
{
    return network[index].get_int_angle(xp,yp);
}

double* actin_ensemble::get_position(int index)
{
    return network[index].getposcm();
}

double actin_ensemble::get_alength(int index)
{
    return network[index].get_length();
}

double* actin_ensemble::get_ends(int index)
{
    return network[index].getendpts();
}

void actin_ensemble::update()
{
    ///Maybe change 6 to 4 for 2d
    av_vel=0;
    for (unsigned int i=0; i<network.size(); i++) {
        vpar=(network[i].get_forces()[0])/network[i].get_friction()[0]  + sqrt(6*temperature/(dt*network[i].get_friction()[0]))*rng_n(0,1);
        vperp=(network[i].get_forces()[1])/network[i].get_friction()[1] + sqrt(6*temperature/(dt*network[i].get_friction()[1]))*rng_n(0,1);
        vx=vpar*cos(network[i].getpos()[2])-vperp*sin(network[i].getpos()[2]);
        vy=vpar*sin(network[i].getpos()[2])+vperp*cos(network[i].getpos()[2]);
        omega=network[i].get_forces()[2]/network[i].get_friction()[2] + sqrt(6*temperature/(dt*network[i].get_friction()[2]))*rng_n(0,1);

        alength=network[i].get_length();

        xnew=network[i].getposcm()[0]+dt*vx;
        ynew=network[i].getposcm()[1]+dt*vy;
        phinew=network[i].getpos()[2]+dt*omega;


        a_ends[0]=xnew-alength*0.5*cos(phinew);
        a_ends[1]=ynew-alength*0.5*sin(phinew);
        a_ends[2]=xnew+alength*0.5*cos(phinew);
        a_ends[3]=ynew+alength*0.5*sin(phinew);

        if (a_ends[0]<=-fov[0]*0.5 || a_ends[0]>=fov[0]*0.5 || a_ends[2]<=-fov[0]*0.5 || a_ends[2]>=fov[0]*0.5)
        {
            vx=-vx;//xnew=network[i].getposcm()[0]-dt*vx;//  
            omega=-omega;//phinew=network[i].getpos()[2]-dt*omega;//omega=-omega;

        }
        if (a_ends[1]<=-fov[1]*0.5 || a_ends[1]>=fov[1]*0.5 || a_ends[3]<=-fov[1]*0.5 || a_ends[3]>=fov[1]*0.5)
        {
            vy=-vy;//ynew=network[i].getposcm()[1]-dt*vy;
            omega=-omega;//phinew=network[i].getpos()[2]-dt*omega;
        }

        xnew=network[i].getposcm()[0]+dt*vx;
        ynew=network[i].getposcm()[1]+dt*vy;
        phinew=network[i].getpos()[2]+dt*omega;
        network.at(i)=actin(xnew,ynew,phinew,alength,fov[0],fov[1],nq[0],nq[1],visc);
    }


}

void actin_ensemble::update_forces(int index, double f1, double f2, double f3)
{
    network[index].update_force(f1,f2,f3);
}


void actin_ensemble::write(std::ofstream& fout)
{
    for (unsigned int i=0; i<network.size(); i++) {
        fout<<network.at(i).getendpts()[0]<<"\t"<<network.at(i).getendpts()[1]<<"\t"<<network.at(i).getendpts()[2]-network.at(i).getendpts()[0]<<"\t"<<network.at(i).getendpts()[3]-network.at(i).getendpts()[1]<<"\n";
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
    double sum = 0, phi1, phi2;
    std::vector<double> corr;
    std::vector<int> monos = mono_map[polymer_index]; 
    for(unsigned int i = 0; i < monos.size() - 1; i++){

        phi1 = network.at(monos[i]).getposcm()[2]; 
        phi2 = network.at(monos[i+1]).getposcm()[2]; 

        sum += cos(phi2 - phi1);

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
    double L = s * ld, sum = 0, phi, sk;

    for (int i = 0; i < s; i++){

        phi = network[ mons[ i ] ].getposcm()[2];
        sk = i * ld + ld/2;
        sum += phi * ld * cos( n * pi * sk / L);

    }

    return sqrt(2.0/L)*sum;

}

// Adds Links between monomers
void actin_ensemble::actin_ensemble::connect_polymers(link_ensemble * links, double link_length, double stretching_stiffness, 
        double bending_stiffness, std::string link_color){
    
    int monomer_index;
    for (unsigned int i = 0; i < mono_map.size(); i++){
        for (unsigned int j = 0; j < mono_map[i].size(); j++){
            monomer_index = mono_map[i][j];
            Link * l = new Link( link_length, stretching_stiffness, bending_stiffness, this, monomer_index, monomer_index + 1,  link_color); 
            links->add_link( l );
            actin_link_map[monomer_index][monomer_index + 1] = l;

        }
    }
}

// Update bending forces between monomers
void actin_ensemble::bending_update(link_ensemble * links){
   
    Link * lft_lnk, ctr_lnk, rt_lnk;
    double forcex, forcey, force_par, force_perp, left_torque, right_torque;
    std::vector<std::vector<double> > forces_x, forces_y;
    

    //initialize all forces to be 0
    for (unsigned int i = 0; i < mono_map.size(); i++){
        for (unsigned int j = 0; j < mono_map[i].size(); j++){
            force_x[i][j] = 0;
            force_y[i][j] = 0;
        }
    }
    
    //Calculate the force at each Link position as outlined by Nedelec, Foethke (2007)
    //Keep the forces at the ends of each filament 0
    for (unsigned int i = 0; i < mono_map.size(); i++){
        
        for (unsigned int j = 1; j < mono_map[i].size() - 1; j++){
            
             
            lft_lnk = actin_link_map[mono_map[i][j-1]][mono_map[i][j  ]];
            ctr_lnk = actin_link_map[mono_map[i][j  ]][mono_map[i][j+1]];
            rt_lnk  = actin_link_map[mono_map[i][j+1]][mono_map[i][j+2]];
            
            forcex = lft_link->get_kb * lft_lnk->get_posx - 2 * ctr_lnk->get_kb * ctr_lnk->get_posx + rt_lnk->get_kb * rt_lnk->get_posx;
            forcey = lft_link->get_kb * lft_lnk->get_posy - 2 * ctr_lnk->get_kb * ctr_lnk->get_posy + rt_lnk->get_kb * rt_lnk->get_posy;

            forces_x[i][j-1] += -1 * force_x;
            forces_x[i][j]   +=  2 * force_x;
            forces_x[i][j+1] += -1 * force_x;

            forces_y[i][j-1] += -1 * force_y;
            forces_y[i][j]   +=  2 * force_y;
            forces_y[i][j+1] += -1 * force_y;
           
            lft_lnk = ctr_lnk;
            ctr_lnk = rt_lnk;
        }
    }

    //Calculate the forces at each center of mass of the monomers
    //Update the monomer
    for (unsigned int i = 0; i < actin_network->mono_map.size(); i++){
        for (unsigned int j = 0; j < actin_network->mono_map[i].size(); j++){
            
            forcex = (forces_x[i][j] + forces_x[i][j+1]) / 2;
            forcey = (forces_y[i][j] + forces_y[i][j+1]) / 2;
            
            force_par   =    forcex*this->get_direction(mono_map[i][j])[0] + forcey*this->get_direction(mono_map[i][j])[1];
            force_perp  =   -forcex*this->get_direction(mono_map[i][j])[1] + forcey*this->get_direction(mono_map[i][j])[0];
            
            if (j == 0)
                left_torque = 0;
            else
                left_torque = cross(left_linkx-this->get_position(mono_map[i][j])[0], left_linky-this->get_position(mono_map[i][j])[1], forcex, forcey);
            
            if (j == mono_map[i])
                right_torque = 0;
            else
                right_torque = cross(right_linkx-this->get_position(mono_map[i][j])[0], right_linky-this->get_position(mono_map[i][j])[1], forcex, forcey);

            actin_network->update_forces(mono_map[i][j], force_par, force_perp, left_torque + right_torque);
            

        }
    }

}


