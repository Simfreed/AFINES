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
    std::vector<double> link_map_value; // xcm, ycm, angle
    double * rod_end_pts;
    double  xcm, ycm, theta;
    for (int i=0; i<npolymer; i++) {
        theta=3*pi/2;//rng(0,2*pi);
        //nmonomer = (int) rng(nmonomer_min, nmonomer_max);
        //the start of the polymer: 
        //            network.push_back(actin(rng(-0.5*(view*fovx-ld),0.5*(view*fovx-ld)), rng(-0.5*(view*fovy-ld),0.5*(view*fovy-ld)),
        //                        theta,ld,fov[0],fov[1],nq[0],nq[1],visc));
        std::cout<<"WARNING: STARTING ACTIN FILAMENT POSITION CHOSEN DETERMINISTICALLY\n";
        network.push_back(actin(0,0,theta,ld,fov[0],fov[1],nq[0],nq[1],visc));
        //Add the quadrants of the first rod
        std::vector<std::vector<int> > tmp_quads=network.back().get_quadrants();
        for (int xindex=0; xindex<tmp_quads[0].size(); xindex++) {
            for (int yindex=0; yindex<tmp_quads[1].size(); yindex++) {
                quad_fils[tmp_quads[0][xindex]][tmp_quads[1][yindex]].push_back(network.size()-1);
            }
        }
        // add monomers to the polymer
        mono_map[i] = empty_vector; 
        for (int j=0; j<nmonomer-1; j++) {

            // Calculate a link (represented as a dead motor) at the end of the rod
            ltheta = theta; //rng(0, 2*pi);
            rod_end_pts = network.back().getendpts();

            lx=rod_end_pts[2] + 0.5*link_ld*cos(ltheta);
            ly=rod_end_pts[3] + 0.5*link_ld*sin(ltheta);

            // Calculate the Next rod on the actin polymer--  continues from the static motor
            //theta = rng(0,2*pi);
            xcm = lx + link_ld*0.5*cos(ltheta) + ld*0.5*cos(theta);
            ycm = ly + link_ld*0.5*sin(ltheta) + ld*0.5*sin(theta);

            // Check that this monomer is in the field of view, otherwise start a new polymer:
            if ( xcm > (0.5*(view*fovx - ld)) || xcm < (-0.5*(view*fovx - ld)) 
                    || ycm > (0.5*(view*fovy - ld)) || ycm < (-0.5*(view*fovy - ld)) )
            {
                std::cout<<"\nDEBUG:"<<j+1<<"th monomer of "<<i<<"th polymer outside field of view; stopped building polymer\n";
                break;
            }else{

                // Load the link map: {monomer index} --> {link_xcm, link_ycm, link_angle} 
                link_map_value.clear();
                link_map_value.push_back(lx);
                link_map_value.push_back(ly);
                link_map_value.push_back(ltheta);
                link_map[ network.size()-1 ] = link_map_value;
                mono_map[ i ].push_back(network.size() -1);

                // Add the actin monomer
                network.push_back( actin(xcm, ycm, theta, ld, fov[0], fov[1], nq[0], nq[1], visc) );

                // Add it's quadrants:
                std::vector<std::vector<int> > tmp_quads=network.back().get_quadrants();
                for (int xindex=0; xindex<tmp_quads[0].size(); xindex++) {
                    for (int yindex=0; yindex<tmp_quads[1].size(); yindex++) {
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
    for (int i=0; i<network.size(); i++) {
        std::vector<std::vector<int> > tmp_quads=network[i].get_quadrants();
        for (int xindex=0; xindex<tmp_quads[0].size(); xindex++) {
            for (int yindex=0; yindex<tmp_quads[1].size(); yindex++) {
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
    av_vel=0;
    for (int i=0; i<network.size(); i++) {
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
    for (int i=0; i<network.size(); i++) {
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
    for(int i = 0; i < monos.size() - 1; i++){

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

    for(int i = 0; i < mono_map.size(); i++){
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

// Adds dead motors between monomers of an actin polymer to act as springs.
// Needs a motor_ensemble object to work, first 
void actin_ensemble::actin_ensemble::connect_polymers(link_ensemble * links, double link_length, double link_stiffness, 
        std::string link_color){
    int monomer_index;
    std::vector<double> motor_coords;
    for (int i = 0; i < mono_map.size(); i++){
        for (int j = 0; j < mono_map[i].size(); j++){
            monomer_index = mono_map[i][j];
            //            motor_coords = link_map[monomer_index];
            links->add_link( Link( link_length, link_stiffness, this, monomer_index, monomer_index + 1,  link_color) );
        }
    }
}
// Adds dead motors between monomers of an actin polymer to act as springs.
// Add additional dead motors between center of masses of monomers to account for bending energy
// Needs a motor_ensemble object to work, first 
void actin_ensemble::actin_ensemble::connect_polymers(link_ensemble * links, 
        double link_length, double link_stiffness, std::string link_color,
        double b_link_stiffness, std::string b_link_color){

    int mono1, mono2;
    double x1, x2, y1, y2, ang, link_length2, b_link_length1, b_link_length2;
    std::vector<double> motor_coords;
    for (int i = 0; i < mono_map.size(); i++){
        for (int j = 0; j < mono_map[i].size(); j++){

            // Join monomers within the polymer
            mono1 = mono_map[i][j];
            mono2 = mono1 + 1; //mono_map[i][j+1];

            link_length2 = 10; //dis_points(x1, x2, y1, y2) + 100 * eps;
            std::cout<<"DEBUG: motor angle = "<<ang<<"\n";
            links->add_link( Link( link_length, link_stiffness, this, mono1, mono2,  link_color) );

            // Add bending energy links
            x1      = network[mono1].getposcm()[0];
            y1      = network[mono1].getposcm()[1];
            x2      = network[mono2].getposcm()[0];
            y2      = network[mono2].getposcm()[1];
            b_link_length1 = network[mono1].get_length()/2 + network[mono2].get_length()/2 + link_length;
            b_link_length2 = dis_points(x1,y1,x2,y2);
            //            std::cout<<"DEBUG: b_link_length2-b_link_length1 : "<<b_link_length2-b_link_length1<<"\n";

            //            mots->add_motor( motor( (x2+x1)/2, (y2+y1)/2, motor_coords[2],//atan2(y2-y1,x2-x1), 
            //                        b_link_length1, this, 1, 1, mono1, mono2, fov[0], fov[1], 0, b_link_stiffness,
            //                        0, 0, 0, ld, visc, b_link_color) ); 

        }
    }
}
