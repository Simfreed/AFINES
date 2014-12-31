/*
 *  filament.cpp
 *  
 *
 *  Created by Simon Freedman on 12/22/2014
 *  Copyright 2014 University of Chicago. All rights reserved.
 *
 */

#include "filament.h"
#include "actin.h"
#include "globals.h"

filament::filament(){}

filament::filament(double startx, double starty, double startphi, int nrod, double fovx, double fovy, int nqx, int nqy, 
        double visc, double deltat, double temp, bool isStraight,
        double rodLength, double linkLength, double stretching_stiffness, double bending_stiffness,
        double frac_force)
{
    fov[0] = fovx;
    fov[1] = fovy;
    nq[0] = nqx;
    nq[1] = nqy;
    dt = deltat;
    temperature = temp;
    gamma = 0;
    fracture_force = frac_force;

    //the start of the polymer: 
    rods.push_back( new actin( startx, starty, startphi, rodLength, fov[0], fov[1], nq[0], nq[1], visc) );
    lks.push_back( new Link(linkLength, stretching_stiffness, bending_stiffness, this, -1, 0) );  
    
    double  xcm, ycm, lphi, phi;
    phi = startphi;
    for (int j = 1; j < nrod; j++) {

        // Calculate the Next rod on the actin polymer--  continues from the link
        if (!isStraight){ 

            //constrain angle between consecutive segments to be small because that's where
            //Nedelec/ Foethke claim their bending energy regime matters
            
            phi += rng(-1*maxSmallAngle , maxSmallAngle);
        }
            
        lphi = (phi + rods.back()->get_angle())/2;
        xcm = rods.back()->get_end()[0] + linkLength*cos(lphi) + rodLength*0.5*cos(phi);
        ycm = rods.back()->get_end()[1] + linkLength*sin(lphi) + rodLength*0.5*sin(phi);

        // Check that this monomer is in the field of view; if not stop building the polymer
        if (       xcm > (0.5*(fov[0] - rodLength)) || xcm < (-0.5*(fov[0] - rodLength)) 
                || ycm > (0.5*(fov[1] - rodLength)) || ycm < (-0.5*(fov[1] - rodLength))      )
        {
            std::cout<<"DEBUG:"<<j+1<<"th segment of filament outside field of view; stopped building filament\n";
            break;
        }else{
            // Add the segment
            rods.push_back( new actin(xcm, ycm, phi, rodLength, fov[0], fov[1], nq[0], nq[1], visc) );
            lks.push_back( new Link(linkLength, stretching_stiffness, bending_stiffness, this, j-1, j) );  
        } 

    }
   
    lks.push_back( new Link(linkLength, stretching_stiffness, bending_stiffness, this, rods.size()-1, -1) );  
    tostring = this->to_string();
}

filament::filament(std::vector<actin *> rodvec, double linkLength, double stretching_stiffness, double bending_stiffness, 
        double deltat, double temp, double frac_force, double g){

    fov[0] = rodvec[0]->get_fov()[0];
    fov[1] = rodvec[0]->get_fov()[1];
    nq[0] = rodvec[0]->get_nq()[0];
    nq[1] = rodvec[0]->get_nq()[1];
    dt = deltat;
    //rods = rodvec;
    temperature = temp;
    fracture_force = frac_force;
    gamma = g;
    
    //Link em up
    for (unsigned int j = 0; j < rodvec.size(); j++) {

        rods.push_back(new actin(*(rodvec[j])));
        lks.push_back( new Link(linkLength, stretching_stiffness, bending_stiffness, this, j-1, j) );  

    }

    if (rods.size() > 0){
        lks.push_back( new Link(linkLength, stretching_stiffness, bending_stiffness, this, rods.size() - 1, -1) );  
    }
    tostring = this->to_string();
}

filament::~filament(){
    
    int nr = rods.size(), nl = lks.size();
    for (int i = 0; i < nr; i ++)
        delete rods[i];
    for (int i = 0; i < nl; i ++)
        delete lks[i];
    
    rods.clear();
    lks.clear();
}

std::vector<std::vector<std::vector<int> > > filament::get_quadrants()
{
    //should return a map between rod and x, y coords of quadrant
    std::vector<std::vector<std::vector<int> > > quads;
    
    for (unsigned int i=0; i<rods.size(); i++) {
        
        quads.push_back(rods[i]->get_quadrants());
    }
    return quads;
}


void filament::update(double t)
{
    double vpar, vperp, vx, vy, omega, alength, xnew, ynew, phinew, phiprev, a_ends[4]; 
    double xleft, xright;
   
    for (unsigned int i = 0; i < rods.size(); i++){

        double * fric = rods[i]->get_friction();
        vpar=(rods[i]->get_forces()[0])/fric[0]  + sqrt(2*temperature/(dt*fric[0]))*rng_n(0,1);
        vperp=(rods[i]->get_forces()[1])/fric[1] + sqrt(2*temperature/(dt*fric[1]))*rng_n(0,1);
        vx=vpar*cos(rods[i]->get_angle())-vperp*sin(rods[i]->get_angle());
        vy=vpar*sin(rods[i]->get_angle())+vperp*cos(rods[i]->get_angle());
        omega=rods[i]->get_forces()[2]/fric[2] + sqrt(2*temperature/(dt*fric[2]))*rng_n(0,1);
        delete[] fric;

        alength=rods[i]->get_length();

        xnew=rods[i]->get_xcm()+dt*vx;
        ynew=rods[i]->get_ycm()+dt*vy;
        phinew=rods[i]->get_angle()+dt*omega;


        a_ends[0]=xnew-alength*0.5*cos(phinew);
        a_ends[1]=ynew-alength*0.5*sin(phinew);
        a_ends[2]=xnew+alength*0.5*cos(phinew);
        a_ends[3]=ynew+alength*0.5*sin(phinew);

        //Calculate the sheared simulation bounds (at this height)
        xleft  = std::max(-fov[0] * 0.5 + gamma * a_ends[1] * t, -fov[0] * 0.5 + gamma * a_ends[3] * t);
        xright = std::min( fov[0] * 0.5 + gamma * a_ends[1] * t,  fov[0] * 0.5 + gamma * a_ends[3] * t);

        if (a_ends[0]<=xleft || a_ends[0]>=xright || a_ends[2]<=xleft || a_ends[2]>=xright)
        {
            vx=-vx;//xnew=rods[i]->get_xcm()-dt*vx;//  
            omega=-omega;//phinew=rods[i]->get_angle()-dt*omega;//omega=-omega;

        }
        if (a_ends[1]<=-fov[1]*0.5 || a_ends[1]>=fov[1]*0.5 || a_ends[3]<=-fov[1]*0.5 || a_ends[3]>=fov[1]*0.5)
        {
            vy=-vy;//ynew=rods[i]->get_ycm()-dt*vy;
            omega=-omega;//phinew=rods[i]->get_angle()-dt*omega;
        }

        xnew=rods[i]->get_xcm()+dt*vx;
        ynew=rods[i]->get_ycm()+dt*vy;
        phinew=rods[i]->get_angle()+dt*omega;

        // Keep consecutive angles small 

        if (i >= 1){

            phiprev = rods[i-1]->get_angle();

            if( phinew - phiprev > maxSmallAngle )
            {
                phinew = phiprev + maxSmallAngle;
            }
            if( phinew - phiprev < -1 * maxSmallAngle)
            {
                phinew = phiprev - maxSmallAngle;
            }

        }
        rods[i]->set_xcm(xnew);
        rods[i]->set_ycm(ynew);
        rods[i]->set_phi(phinew);
        rods[i]->update(); //updates all derived quantities (e.g., endpoints, forces = 0, etc.)

    }
}

void filament::update_bending()
{

    double forcex, forcey, force_par, force_perp, lft_trq, rt_trq;
    std::vector<double> node_forces_x, node_forces_y;
    
    if (rods.size() < 2){ //no bending energy!
        return;
    }
    
    //initialize all NODE forces to be 0
    for (unsigned int j = 0; j <= rods.size(); j++){
        node_forces_x.push_back(0);
        node_forces_y.push_back(0);
    }
    
    //Calculate the force at each Link position as outlined by Nedelec, Foethke (2007)
    //Keep the forces at the ends of each filament 0
    
    for (unsigned int j = 2; j <= rods.size(); j++){

        forcex =    lks[j-2]->get_kb() * lks[j-2]->get_posx() 
              - 2 * lks[j-1]->get_kb() * lks[j-1]->get_posx() 
              +     lks[ j ]->get_kb() * lks[ j ]->get_posx(); 
        
        forcey =    lks[j-2]->get_kb() * lks[j-2]->get_posy() 
              - 2 * lks[j-1]->get_kb() * lks[j-1]->get_posy() 
              +     lks[ j ]->get_kb() * lks[ j ]->get_posy();

        node_forces_x[j-2] += -1 * forcex;
        node_forces_x[j-1] +=  2 * forcex;
        node_forces_x[j  ] += -1 * forcex;

        node_forces_y[j-2] += -1 * forcey;
        node_forces_y[j-1] +=  2 * forcey;
        node_forces_y[j  ] += -1 * forcey;
    
    }
    
    //Calculate the forces at each center of mass of the segments
    //Update the monomer

    for (unsigned int j = 0; j < rods.size(); j++){
    
        forcex = (node_forces_x[j] + node_forces_x[j+1]) / 2;
        forcey = (node_forces_y[j] + node_forces_y[j+1]) / 2;

        force_par   =  forcex*rods[j]->get_direction()[0] + forcey*rods[j]->get_direction()[1];
        force_perp  = -forcex*rods[j]->get_direction()[1] + forcey*rods[j]->get_direction()[0];
        
        if (j == 0)
            lft_trq = 0;
        else
            lft_trq = cross(rods[j]->get_xcm() - lks[j]->get_posx(),
                            rods[j]->get_ycm() - lks[j]->get_posy(), forcex, forcey);

        if (j == rods.size() - 1)
            rt_trq = 0;
        else{
            rt_trq = cross(rods[j]->get_xcm() - lks[j+1]->get_posx(),
                           rods[j]->get_ycm() - lks[j+1]->get_posy(), forcex, forcey);
        }
        
        rods[j]->update_force(force_par, force_perp, lft_trq + rt_trq);

    }

}

std::vector<filament *> filament::update_stretching()
{
    std::vector<filament *> newfilaments;
    for (unsigned int i=0; i < lks.size(); i++) {
        lks[i]->step();
        if (fabs(lks[i]->get_stretch_force()) > fracture_force){
            newfilaments = this->fracture(i);
            break;
        }
        else
            lks[i]->filament_update();
    }
    return newfilaments;
}

actin * filament::get_rod(int i)
{
    return rods[i];
}

Link * filament::get_link(int i)
{
    return lks[i];
}

void filament::update_shear(){
    
    double ycm, forcex, force_par, force_perp, lft_trq, rt_trq;
    
    for (unsigned int i = 0; i < rods.size(); i++){
        
        ycm = rods[i]->get_ycm();
        forcex = gamma * ycm;
        
        if (i == 0)
            lft_trq = 0;
        else
            lft_trq = -1 * forcex * (ycm - lks[i]->get_posy()); //cross product with fy = 0

        if (i == rods.size() - 1)
            rt_trq = 0;
        else{
            rt_trq = -1 * forcex * (ycm - lks[i+1]->get_posy());
        }
        
        force_par   =  forcex*rods[i]->get_direction()[0];
        force_perp  = -forcex*rods[i]->get_direction()[1];
        
        rods[i]->update_force(force_par, force_perp, lft_trq + rt_trq);
    }


}

void filament::update_forces(int index, double f1, double f2, double f3)
{
    rods[index]->update_force(f1,f2,f3);
}

void filament::set_shear(double g){

    gamma = g;

}

std::string filament::write_rods(){
    std::string all_rods;
    for (unsigned int i =0; i < rods.size(); i++)
    {
        all_rods += rods[i]->write();
    }

    return all_rods;
}

std::string filament::write_links(){
    std::string all_links;
    for (unsigned int i =0; i < lks.size(); i++)
    {
        all_links += lks[i]->write();
    }

    return all_links;
}

std::vector<actin *> filament::get_rods(unsigned int first, unsigned int last)
{
    std::vector<actin *> newrods;
    for (unsigned int i = first; i <= last; i++)
    {
        if (i >= rods.size())
            break;
        else
            newrods.push_back(rods[i]);
    }
    return newrods;
}

std::vector<filament *> filament::fracture(int node){

    std::vector<filament *> newfilaments;
    std::cout<<"DEBUG: fracturing";

    newfilaments.push_back(
            new filament(this->get_rods(0,           node - 1), lks[0]->get_length(), lks[0]->get_kl(), lks[0]->get_kb(), 
                dt, temperature, fracture_force, gamma));
    newfilaments.push_back(
            new filament(this->get_rods(node, rods.size() - 1), lks[0]->get_length(), lks[0]->get_kl(), lks[0]->get_kb(), 
                dt, temperature, fracture_force, gamma));

    return newfilaments;

}

bool filament::operator==(const filament& that){
    
    for (unsigned int i = 0; i < rods.size(); i++)
        if (!(*(rods[i]) == *(that.rods[i])))
            return false;
    
    for (unsigned int i = 0; i < lks.size(); i++)
        if (!(lks[i]->is_similar(*(that.lks[i]))))
            return false;

    return (this->fov[0] == that.fov[0] && this->fov[1] == that.fov[1] && 
            this->nq[0] == that.nq[0] && this->nq[1] == that.nq[1] &&
            this->gamma == that.gamma && this->temperature == that.temperature &&
            this->dt == that.dt && this->fracture_force == that.fracture_force);

}

std::string filament::to_string(){
    
    // Note: not including links in to_string, because link's to_string includes filament's to_string
    char buffer[200];
    std::string out = "";

    for (unsigned int i = 0; i < rods.size(); i++)
        out += rods[i]->to_string();

    sprintf(buffer, "fov = (%f, %f)\tnq = (%d, %d)\tgamma = %f\ttemperature = %f\tdt = %f\tfracture_force=%f\n",
            fov[0], fov[1], nq[0], nq[1], gamma, temperature, dt, fracture_force);
   
    return out + buffer; 

}
