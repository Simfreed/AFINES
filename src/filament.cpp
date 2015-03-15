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
//using namespace std;
filament::filament(){

    fov[0] = 50;
    fov[1] = 50;
    nq[0] = 100;
    nq[1] = 100;
    dt = 0.001;
    temperature = 0;
    gamma = 0;
    fracture_force = 1000000;
    BC = "REFLECTIVE";

}

filament::filament(double startx, double starty, double startphi, int nactin, double fovx, double fovy, int nqx, int nqy, 
        double visc, double deltat, double temp, bool isStraight,
        double ballRadius, double linkLength, double stretching_stiffness, double bending_stiffness,
        double frac_force, string bdcnd)
{
 
    fov[0] = fovx;
    fov[1] = fovy;
    nq[0] = nqx;
    nq[1] = nqy;
    dt = deltat;
    temperature = temp;
    gamma = 0;
    fracture_force = frac_force;
    BC = bdcnd;
    kb = bending_stiffness;

    //the start of the polymer: 
    actins.push_back(new actin( startx, starty, ballRadius, fov[0], fov[1], nq[0], nq[1], viscosity));
    
    double  xcm, ycm, phi, variance;
    phi = startphi;
    
    if (temp != 0) variance = temp/(bending_stiffness * linkLength * linkLength);
    else variance = 1;

    for (int j = 1; j < nactin; j++) {

        xcm = actins.back()->get_xcm() + linkLength*cos(phi);
        ycm = actins.back()->get_ycm() + linkLength*sin(phi);

        // Check that this monomer is in the field of view; if not stop building the polymer
        if (       xcm > (0.5*(fov[0] - ballRadius)) || xcm < (-0.5*(fov[0] - ballRadius)) 
                || ycm > (0.5*(fov[1] - ballRadius)) || ycm < (-0.5*(fov[1] - ballRadius))      )
        {
            cout<<"DEBUG:"<<j+1<<"th segment of filament outside field of view; stopped building filament\n";
            break;
        }else{
            // Add the segment
            actins.push_back( new actin(xcm, ycm, phi, ballRadius, fov[0], fov[1], nq[0], nq[1], viscosity) );
            links.push_back( new Link(linkLength, stretching_stiffness, this, j-1, j) );  
            
        } 
        
        // Calculate the Next angle on the actin polymer
        if (!isStraight){ 

            //phi += rng(-1*maxSmallAngle , maxSmallAngle);
            phi += rng_n(0, variance);
        }

    }
   
}

filament::filament(vector<actin *> actinvec, double linkLength, double stretching_stiffness, double bending_stiffness, 
        double deltat, double temp, double frac_force, double g, string bdcnd)
{
    
    dt = deltat;
    temperature = temp;
    fracture_force = frac_force;
    gamma = g; 
    BC = bdcnd;
    kb = bending_stiffness;

    if (actinvec.size() > 0)
    {
        fov[0] = actinvec[0]->get_fov()[0];
        fov[1] = actinvec[0]->get_fov()[1];
        nq[0] = actinvec[0]->get_nq()[0];
        nq[1] = actinvec[0]->get_nq()[1];

        actins.push_back(new actin(*(actinvec[0])));
    }
    else{
        fov[0] = 0;
        fov[1] = 0;
        nq[0] = 0;
        nq[1] = 0;
    }

    //Link em up
    for (unsigned int j = 1; j < actinvec.size(); j++) {
        
        links.push_back( new Link(linkLength, stretching_stiffness, bending_stiffness, this, j-1, j) );  
        actins.push_back(new actin(*(actinvec[j])));
    
    }

}

filament::~filament(){
    
    int nr = actins.size(), nl = links.size();
    for (int i = 0; i < nr; i ++)
    {    
        //cout<<"\nDEBUG: deleting pointer "<<actins[i];     
        delete actins[i];
    }
    for (int i = 0; i < nl; i ++)
        delete links[i];
    
    actins.clear();
    links.clear();
}

void filament::add_actin(actin * a, double linkLength, double stretching_stiffness){
    
    actins.push_back(new actin(*a));
    
    if (actins.size() > 1)
        links.push_back( new Link(linkLength, stretching_stiffness, this, actins.size()-1,  actins.size() ) );  
    
}

vector<vector<array<int,2> > > > filament::get_quadrants()
{
    //should return a map between actin and x, y coords of quadrant
    vector<vector<array<int,2> > > quads;
    
    for (unsigned int i=0; i<links.size(); i++) 
        quads.push_back(links[i]->get_quadrants());
    
    return quads;
}

void filament::update(double t)
{
    double vpar, vperp, vx, vy, omega, alength, xnew, ynew, phinew, a_ends[4]; 
    double xleft = -fov[0]*0.5, xright=fov[0]*0.5;
    bool recenter_filament = false;
    double xleft  = -fov[0] * 0.5;
    double xright =  fov[0] * 0.5;
    double yleft  = -fov[1] * 0.5;
    double yright =  fov[1] * 0.5;

    for (unsigned int i = 0; i < actins.size(); i++){
        
        double fric = actins[i]->get_friction();
        vx  = (actins[i]->get_forces()[0])/fric  + sqrt(2*temperature/(dt*fric))*rng_n(0,1);
        vy  = (actins[i]->get_forces()[1])/fric  + sqrt(2*temperature/(dt*fric))*rng_n(0,1);

        xnew = actins[i]->get_xcm()+dt*vx;
        ynew = actins[i]->get_ycm()+dt*vy;

        //Calculate the sheared simulation bounds (at this height)
        if (gamma != 0){
            xleft  = -fov[0] * 0.5 + gamma * ynew[1] * t;
            xright =  fov[0] * 0.5 + gamma * ynew[1] * t;
        }

        if(this->get_BC() == "REFLECTIVE")
        {
            if (xnew <= xleft || xnew >= xright) vx=-vx;
            if (ynew <= yleft || ynew >= yright) vy=-vy;
            
            xnew=actins[i]->get_xcm()+dt*vx;
            ynew=actins[i]->get_ycm()+dt*vy;
        
        }
        else if(this->get_BC() == "PERIODIC")
        {
            if (xnew < xleft) xnew += fov[0];
            else if (xnew > xright) xnew -= fov[0];
         
            
            if (ynew < yleft) ynew += fov[1];
            else if(ynew > yright) ynew -= fov[1];

        }
        
        else if(this->get_BC() == "NONE")
        {
            recenter_filament = true;
        }
        
        actins[i]->set_xcm(xnew);
        actins[i]->set_ycm(ynew);
        actins[i]->update(); //updates all derived quantities (e.g., endpoints, forces = 0, etc.)
    }
    
    if (recenter_filament && actins.size() > 0){
        double midx = actins[floor(actins.size()/2)]->get_xcm();
        double midy = actins[floor(actins.size()/2)]->get_ycm();
        for (unsigned int i = 0; i < actins.size(); i++){
            
            xnew = actins[i]->get_xcm()-midx;
            ynew = actins[i]->get_ycm()-midy;
            
            actins[i]->set_xcm(xnew);
            actins[i]->set_ycm(ynew);
            actins[i]->update();
       
        }
    }

}

vector<filament *> filament::update_stretching()
{
    vector<filament *> newfilaments;
    
    if(links.size() == 0)
        return newfilaments;
    
    for (unsigned int i=0; i < links.size(); i++) {
        links[i]->step();
        if (fabs(links[i]->get_stretch_force()) > fracture_force){
            newfilaments = this->fracture(i);
            break;
        }
        else 
            links[i]->filament_update;
    }
    
    return newfilaments;
}


actin * filament::get_actin(int i)
{
    try
    {
        return actins[i];
    }
    catch (int e)
    {
        cout<<"\nDEBUG: an exception occured while returning the actins[ "<<i<<"]";
        actin * a;
        return a;
    }
}

Link * filament::get_link(int i)
{
    return links[i];
}

void filament::update_shear(){
    
    for (unsigned int i = 0; i < actins.size(); i++)
        actins[i]->update_force( gamma * actins[i]->get_ycm() , 0);
    
}

void filament::update_forces(int index, double f1, double f2)
{
    actins[index]->update_force(f1,f2);
}

void filament::set_shear(double g){
    gamma = g;
}

string filament::write_actins(){
    string all_actins;
    for (unsigned int i =0; i < actins.size(); i++)
    {
        all_actins += actins[i]->write();
    }

    return all_actins;
}

string filament::write_links(){
    string all_links;
    for (unsigned int i =0; i < links.size(); i++)
    {
        all_links += links[i]->write();
    }

    return all_links;
}

vector<actin *> filament::get_actins(unsigned int first, unsigned int last)
{
    vector<actin *> newactins;
    for (unsigned int i = first; i < last; i++)
    {
        if (i < 0 || i >= actins.size())
            break;
        else
            newactins.push_back(new actin(*(actins[i])));
    }
    return newactins;
}

vector<filament *> filament::fracture(int node){

    vector<filament *> newfilaments;
    cout<<"\n\tDEBUG: fracturing at node "<<node;
    
    vector<actin *> lower_half = this->get_actins(0, node);
    vector<actin *> upper_half = this->get_actins(node, actins.size());

    newfilaments.push_back(
            new filament(lower_half, links[0]->get_length(), links[0]->get_kl(), kb, 
                dt, temperature, fracture_force, gamma, BC));
    newfilaments.push_back(
            new filament(upper_half, links[0]->get_length(), links[0]->get_kl(), kb, 
                dt, temperature, fracture_force, gamma, BC));

    int s = actinvec.size();
    for (int i = 0; i < node; i++) delete lower_half[i];
    for (int i = node; i < actins.size(); i++) delete upper_half[i];
    
    lower_half.clear();
    upper_half.clear();
    
    return newfilaments;

}

bool filament::operator==(const filament& that){
    
    for (unsigned int i = 0; i < actins.size(); i++)
        if (!(*(actins[i]) == *(that.actins[i])))
            return false;
    
    for (unsigned int i = 0; i < links.size(); i++)
        if (!(links[i]->is_similar(*(that.links[i]))))
            return false;

    return (this->fov[0] == that.fov[0] && this->fov[1] == that.fov[1] && 
            this->nq[0] == that.nq[0] && this->nq[1] == that.nq[1] &&
            this->gamma == that.gamma && this->temperature == that.temperature &&
            this->dt == that.dt && this->fracture_force == that.fracture_force);

}

string filament::to_string(){
    
    // Note: not including links in to_string, because link's to_string includes filament's to_string
    char buffer[200];
    string out = "";

    for (unsigned int i = 0; i < actins.size(); i++)
        out += actins[i]->to_string();

    sprintf(buffer, "fov = (%f, %f)\tnq = (%d, %d)\tgamma = %f\ttemperature = %f\tdt = %f\tfracture_force=%f\n",
            fov[0], fov[1], nq[0], nq[1], gamma, temperature, dt, fracture_force);
   
    return out + buffer; 

}

string filament::get_BC(){
    return BC; 
}

void filament::set_BC(string s){
    this->BC = s;
}
    
void filament::fwd_bending_update()
{
    //Calculate the force at each Link position as outlined by 
    // Allen and Tildesley in Computer Simulation of Liquids, Appendix C
    // Here we let db = (xb, yb) = (linkx[b] - linkx[b-1], linky[b] - linky[b-1])
 
    double xam1, xa, xap1, xap2;
    double yam1, ya, yap1, yap2;
    double Cam1am1, Caa, Cap1ap1, Cap2ap2, Caam1, Caap1, Cap1ap2;
    double theta_a, theta_ap1, theta_ap2;
    double coef1, coef2, coef3;
    double fx1, fx2, fx3, fy1, fy2, fy3;
    
    double forcex, forcey;
    //initialize all NODE forces to be 0
    
    //single actin --> no bending energy 
    if (actins.size() > 2)
    { 
        //First two actins won't have any bending forces: 
        
        // More than 1 actin--> bending forces to calculate
        xam1      = actins[1]->get_xcm() - actins[0]->get_xcm();
        yam1      = actins[1]->get_ycm() - actins[0]->get_ycm();

        xa        = actins[2]->get_xcm() - actins[1]->get_xcm();
        ya        = actins[2]->get_ycm() - actins[1]->get_ycm();

        Cam1am1 = xam1*xam1 + yam1*yam1;
        Caa     = xa  *xa   + ya  *ya;
        Caam1   = xam1*xa   + yam1*ya;

        theta_a   = links[1]->get_angle() - links[0]->get_angle();

        if( fabs(theta_a) < maxSmallAngle )
            coef1 = -kb * ( 1 / sqrt( Cam1am1 * Caa     ) ); 
        else
            coef1 = -kb*theta_a   / sin(theta_a)   * ( 1 / sqrt( Cam1am1 * Caa     ) ); 

        //2 actins --> bending force on last link; only has one term
        if (actins.size() == 2)
        {

            fx1 = coef1 * ( Caam1/Caa * xa - xam1 );
            fy1 = coef1 * ( Caam1/Caa * ya - yam1 );
            //cout<<"\nDEBUG: magnitude of bending forces at actins: ( "<<fx1<<" , "<<fy1<<" )";
            
            if (fabs(fx1) < eps)
                fx1 = 0;
            if (fabs(fy1) < eps)
                fy1 = 0;

            actins[2]->update_force(fx1, fy1);

        }
        else
        {
            //More than 2 actins--> more bending forces to calculate
            xap1 = actins[3]->get_xcm() - actins[2]->get_xcm();
            yap1 = actins[3]->get_ycm() - actins[2]->get_ycm();

            Cap1ap1 = xap1*xap1 + yap1*yap1;
            Caap1   = xap1*xa   + yap1*ya;

            theta_ap1 = links[2]->get_angle() - links[1]->get_angle();

            if ( fabs(theta_ap1) < maxSmallAngle )
                coef2 =  kb * ( 1 / sqrt( Caa     * Cap1ap1 ) );
            else
                coef2 =  kb*theta_ap1 / sin(theta_ap1) * ( 1 / sqrt( Caa     * Cap1ap1 ) );

            //Enter loop if more than 3 actins. For 3 actin case, the loop is skipped
            for (unsigned int j = 2; j < actins.size() - 2; j++){

                xap2 = actins[j+2]->get_xcm() - actins[j+1]->get_xcm();
                yap2 = actins[j+2]->get_ycm() - actins[j+1]->get_ycm();

                Cap2ap2 = xap2*xap2 + yap2*yap2;
                Cap1ap2 = xap1*xap2 + yap1*yap2;

                theta_ap2 = links[j+1]->get_angle() - links[j]->get_angle();
                
                if ( fabs(theta_ap2) < maxSmallAngle )
                    coef3 = -kb * ( 1 / sqrt( Cap1ap1 * Cap2ap2 ) );
                else
                    coef3 = -kb*theta_ap2 / sin(theta_ap2) * ( 1 / sqrt( Cap1ap1 * Cap2ap2 ) );
                
                fx1 = coef1 * ( Caam1/Caa * xa - xam1 );
                fy1 = coef1 * ( Caam1/Caa * ya - yam1 );

                fx2 = coef2 * ( (1 + Caap1/Cap1ap1) * xap1 - (1 + Caap1/Caa) * xa ); 
                fy2 = coef2 * ( (1 + Caap1/Cap1ap1) * yap1 - (1 + Caap1/Caa) * ya ); 

                fx3 = coef3 * ( xap2 - Cap1ap2/Cap1ap1 * xap1);
                fy3 = coef3 * ( yap2 - Cap1ap2/Cap1ap1 * yap1);

                forcex = fx1 + fx2 + fx3;
                forcey = fy1 + fy2 + fy3;
                if (fabs(forcex) < eps)
                    forcex = 0;
                if (fabs(forcey) < eps)
                    forcey = 0;
                
                actins[i]->update_force(forcex, forcey);

                //increment all variables for next iteration:
                xam1 = xa;
                xa   = xap1;
                xap1 = xap2;

                yam1 = ya;
                ya   = yap1;
                yap1 = yap2;

                Cam1am1 = Caa;
                Caa     = Cap1ap1;
                Cap1ap1 = Cap2ap2;

                Caam1 = Caap1;
                Caap1 = Cap1ap2;

                coef1 = -1*coef2;
                coef2 = -1*coef3;

            }

            //LAST TWO actinS ON THE FILAMENT:
            fx1 = coef1 * ( Caam1/Caa * xa - xam1 );
            fy1 = coef1 * ( Caam1/Caa * ya - yam1 );
            fx2 = coef2 * ( (1 + Caap1/Cap1ap1) * xap1 - (1 + Caap1/Caa) * xa ); 
            fy2 = coef2 * ( (1 + Caap1/Cap1ap1) * yap1 - (1 + Caap1/Caa) * ya ); 
            
            //cout<<"\nDEBUG: 2nd to last actin: (fx1, fy1) + (fx2, fy2) = ( "<<fx1<< " , "<<fy1<<" ) + ( "<<fx2<<" , " <<fy2<<" )";
            
            forcex = fx1 + fx2;
            forcey = fy1 + fy2;
            if (fabs(forcex) < eps)
                forcex = 0;
            if (fabs(forcey) < eps)
                forcey = 0;
            
            actins[ actins.size() - 2 ]->update_force(forcex, forcey);

            /*INCREMENT*/
            xam1 = xa;
            xa   = xap1;

            yam1 = ya;
            ya   = yap1;

            Caa   = Cap1ap1;
            Caam1 = Caap1;

            coef1 = -1*coef2;

            fx1 = coef1 * ( Caam1/Caa * xa - xam1 );
            fy1 = coef1 * ( Caam1/Caa * ya - yam1 );
            //cout<<"\nDEBUG: Last actin: (fx1, fy1) = ( "<<fx1<< " , "<<fy1<<" )";

            if (fabs(fx1) < eps)
                fx1 = 0;
            if (fabs(fy1) < eps)
                fy1 = 0;
            
            actins[ actins.size() - 1 ]->update_force(fx1, fy1);
        }
    }
}


//wrapper, for fwd_bending_update (and bwd bending update if I ever make it)
void filament::update_bending()
{
    if(links.size() > 2)
        this->fwd_bending_update();

}


int filament::get_nactins(){
    return actins.size();
}

double filament::get_bending_energy(){
    
    double sum = 0, theta;

    for (unsigned int i = 0; i < links.size() - 1; i++)
    {
        theta = links[i+1]->get_angle() - links[i]->get_angle();
        sum += theta*theta;
    }
    
    return kb*sum/2.0;

}

double filament::get_stretching_energy(){
    
    if(actins.size() == 0){
        return 0;
    }
    
    double sum = 0, stretch;

    for (unsigned int i = 0; i < links.size(); i++)
    {
        stretch = links[i]->get_stretch_force();
        sum += stretch*stretch;
    }
    
    return sum/(2.0*links[0]->get_l0());

}

double filament::get_total_energy()
{
    return this->get_stretching_energy() + this->get_bending_energy();
}
