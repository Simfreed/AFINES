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
    kinetic_energy = 0;
}

filament::filament(array<double, 2> myfov, array<int, 2> mynq, double deltat, double temp, double shear, 
        double frac, double bending_stiffness, string bndcnd)
{
    fov             = myfov;
    nq              = mynq;
    dt              = deltat;
    temperature     = temp;
    gamma           = shear;
    fracture_force  = frac;
    kb              = bending_stiffness;
    BC              = bndcnd;
    kinetic_energy  = 0;
}

filament::filament(array<double, 3> startpos, int nactin, array<double, 2> myfov, array<int, 2> mynq, double visc, 
        double deltat, double temp, bool isStraight, double actinRadius, double linkLength, double stretching_stiffness,
        double bending_stiffness, double frac_force, string bdcnd)
{
    
    fov = myfov;
    nq = mynq;
    dt = deltat;
    temperature = temp;
    gamma = 0;
    fracture_force = frac_force;
    BC = bdcnd;
    kb = bending_stiffness;
    kinetic_energy  = 0;

    double xcm, ycm, phi, variance;
    //the start of the polymer: 
    actins.push_back(new actin( startpos[0], startpos[1], actinRadius, visc));
    phi = startpos[2];
    
    if (temp != 0) variance = temp/(bending_stiffness * linkLength * linkLength);
    else variance = 1;

    for (int j = 1; j < nactin; j++) {

        xcm = actins.back()->get_xcm() + linkLength*cos(phi);
        ycm = actins.back()->get_ycm() + linkLength*sin(phi);

        // Check that this monomer is in the field of view; if not stop building the polymer
        if (       xcm > (0.5*(fov[0] - actinRadius)) || xcm < (-0.5*(fov[0] - actinRadius)) 
                || ycm > (0.5*(fov[1] - actinRadius)) || ycm < (-0.5*(fov[1] - actinRadius))      )
        {
            cout<<"DEBUG:"<<j+1<<"th segment of filament outside field of view; stopped building filament\n";
            break;
        }else{
            // Add the segment
            actins.push_back( new actin(xcm, ycm, actinRadius, visc) );
            links.push_back( new Link(linkLength, stretching_stiffness, this, {j-1, j}, fov, nq) );  
            
        } 
        
        // Calculate the Next angle on the actin polymer
        if (!isStraight){ 

            //phi += rng(-1*maxSmallAngle , maxSmallAngle);
            phi += rng_n(0, variance);
        }

    }
   
}

filament::filament(vector<actin *> actinvec, array<double, 2> myfov, array<int, 2> mynq, double linkLength, 
        double stretching_stiffness, double bending_stiffness, 
        double deltat, double temp, double frac_force, double g, string bdcnd)
{
    
    dt = deltat;
    temperature = temp;
    fracture_force = frac_force;
    gamma = g; 
    BC = bdcnd;
    kb = bending_stiffness;
    fov = myfov;
    nq = mynq;

    if (actinvec.size() > 0)
    {
        actins.push_back(new actin(*(actinvec[0])));
    }
    
    //Link em up
    if (actinvec.size() > 1){
        for (unsigned int j = 1; j < actinvec.size(); j++) {

            actins.push_back(new actin(*(actinvec[j])));
            links.push_back( new Link(linkLength, stretching_stiffness, this, {(int)j-1, (int)j}, fov, nq) );  

        }
    }
}

filament::~filament(){
    
    cout<<"DELETING FILAMENT\n";
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
    
    if (actins.size() > 1){
        int j = (int) actins.size() - 1;
        links.push_back( new Link(linkLength, stretching_stiffness, this, {j-1,  j}, fov, nq ) );  
    }
}

vector<vector<array<int,2> > > filament::get_quadrants()
{
    //should return a map between actin and x, y coords of quadrant
    vector<vector<array<int,2> > > quads;
    vector<array<int, 2> > l_quads;
    for (unsigned int i=0; i<links.size(); i++){ 
        links[i]->quad_update();
        //cout<<"\nDEBUG links[i]->get_quadrants";
        //l_quads = links[i]->get_quadrants();
        //for_each(l_quads.begin(), l_quads.end(), intarray_printer);
        quads.push_back(links[i]->get_quadrants());
    }
    return quads;
}

void filament::update_positions(double t)
{
    double vx, vy, gamma, T = temperature;
    array<double, 2> newpos;
    kinetic_energy = 0;  

    //if (t < dt*100000) T = 0;
    
    for (unsigned int i = 0; i < actins.size(); i++){
        
        gamma = actins[i]->get_friction();
        
        vx  = (actins[i]->get_force()[0])/gamma  + sqrt(2*T/(dt*gamma))*rng_n(0,1);
        vy  = (actins[i]->get_force()[1])/gamma  + sqrt(2*T/(dt*gamma))*rng_n(0,1);
        
        kinetic_energy += vx*vx + vy*vy;
        
        newpos = boundary_check(i, t, vx, vy); 
        
        actins[i]->set_xcm(newpos[0]);
        actins[i]->set_ycm(newpos[1]);
        actins[i]->reset_force(); 
    }
    
}

array<double, 2> filament::boundary_check(int i, double t, double vx, double vy)
{
    double xnew = actins[i]->get_xcm()+dt*vx;
    double ynew = actins[i]->get_ycm()+dt*vy;
        
    //Calculate the sheared simulation bounds (at this height)
    double xleft  = -fov[0] * 0.5 + gamma * ynew * t; //sheared simulation bounds
    double xright =  fov[0] * 0.5 + gamma * ynew * t;
    double yleft  = -fov[1] * 0.5;
    double yright =  fov[1] * 0.5;

    if(BC == "REFLECTIVE")
    {
        if (xnew <= xleft || xnew >= xright) vx=-vx;
        if (ynew <= yleft || ynew >= yright) vy=-vy;

        xnew=actins[i]->get_xcm()+dt*vx;
        ynew=actins[i]->get_ycm()+dt*vy;

    }
    else if(BC == "PERIODIC")
    {
        if (xnew < xleft)       xnew += fov[0];
        else if (xnew > xright) xnew -= fov[0];


        if (ynew < yleft)      ynew += fov[1];
        else if(ynew > yright) ynew -= fov[1];

    }

    else if(BC == "NONE")
    {
        xnew -= actins[floor(actins.size()/2)]->get_xcm();
        ynew -= actins[floor(actins.size()/2)]->get_ycm();
    }
    
    return {xnew, ynew};

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
            links[i]->filament_update();
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

string filament::write_thermo()
{
    return std::to_string(this->get_kinetic_energy()) + \
        "\t" + std::to_string(this->get_potential_energy()) + \
        "\t" + std::to_string(this->get_total_energy()) + "\n";
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
    
    if(links.size() == 0)
        return newfilaments;

    vector<actin *> lower_half = this->get_actins(0, node+1);
    vector<actin *> upper_half = this->get_actins(node+1, actins.size());

    if (lower_half.size() > 0)
        newfilaments.push_back(
                new filament(lower_half, fov, nq, links[0]->get_length(), links[0]->get_kl(), kb, 
                    dt, temperature, fracture_force, gamma, BC));
    if (upper_half.size() > 0)
        newfilaments.push_back(
                new filament(upper_half, fov, nq, links[0]->get_length(), links[0]->get_kl(), kb, 
                    dt, temperature, fracture_force, gamma, BC));

    for (int i = 0; i < (int)(lower_half.size()); i++) delete lower_half[i];
    for (int i = 0; i < (int)(upper_half.size()); i++) delete upper_half[i];
    
    lower_half.clear();
    upper_half.clear();
    
    return newfilaments;

}

bool filament::operator==(const filament& that){
    
    if (actins.size() != that.actins.size() || links.size() != that.links.size())
        return false;

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

inline double filament::angle_between_links(int i, int j){
    theta = links[i]->get_angle() - links[j]->get_angle();
    return theta - 2*pi*floor(theta/(2*pi)+0.5); // Keep angles between -Pi and Pi
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

        theta_a   = angle_between_links(1,0); //links[1]->get_angle() - links[0]->get_angle();

        if( fabs(theta_a) < maxSmallAngle )
            coef1 = -kb * ( 1 / sqrt( Cam1am1 * Caa     ) ); 
        else
            coef1 = -kb*theta_a   / sin(theta_a)   * ( 1 / sqrt( Cam1am1 * Caa     ) ); 

        //2 actins --> bending force on last link; only has one term
        if (actins.size() == 3)
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

            theta_ap1 = angle_between_links(2,1); //links[2]->get_angle() - links[1]->get_angle();
            
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

                theta_ap2 = angle_between_links(j+1,j); //links[j+1]->get_angle() - links[j]->get_angle();
                
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
                
                actins[j]->update_force(forcex, forcey);

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

void filament::bwd_bending_update()
{
    
    //same thing as fwd_bending_update, but starting from the other end of the filament
    double xam1, xa, xap1, xap2;
    double yam1, ya, yap1, yap2;
    double Cam1am1, Caa, Cap1ap1, Cap2ap2, Caam1, Caap1, Cap1ap2;
    double theta_a, theta_ap1, theta_ap2;
    double coef1, coef2, coef3;
    double fx1, fx2, fx3, fy1, fy2, fy3;
    
    double forcex, forcey;
    //initialize all NODE forces to be 0
   
    int zero = actin.size() - 1, one = actins.size() - 2, two = actins.size()-3, three = actins.size() - 4;
    //single actin --> no bending energy 
    if (actins.size() > 2)
    { 
        //First two actins won't have any bending forces: 
        
        // More than 1 actin--> bending forces to calculate
        xam1      = actins[one]->get_xcm() - actins[zero]->get_xcm();
        yam1      = actins[one]->get_ycm() - actins[zero]->get_ycm();

        xa        = actins[two]->get_xcm() - actins[one]->get_xcm();
        ya        = actins[two]->get_ycm() - actins[one]->get_ycm();

        Cam1am1 = xam1*xam1 + yam1*yam1;
        Caa     = xa  *xa   + ya  *ya;
        Caam1   = xam1*xa   + yam1*ya;

        theta_a   = angle_between_links(one, zero); //links[one]->get_angle() - links[zero]->get_angle();

        if( fabs(theta_a) < maxSmallAngle )
            coef1 = -kb * ( 1 / sqrt( Cam1am1 * Caa     ) ); 
        else
            coef1 = -kb*theta_a   / sin(theta_a)   * ( 1 / sqrt( Cam1am1 * Caa     ) ); 

        //2 actins --> bending force on last link; only has one term
        if (actins.size() == 3)
        {

            fx1 = coef1 * ( Caam1/Caa * xa - xam1 );
            fy1 = coef1 * ( Caam1/Caa * ya - yam1 );
            //cout<<"\nDEBUG: magnitude of bending forces at actins: ( "<<fx1<<" , "<<fy1<<" )";
            
            if (fabs(fx1) < eps)
                fx1 = 0;
            if (fabs(fy1) < eps)
                fy1 = 0;

            actins[two]->update_force(fx1, fy1);

        }
        else
        {
            //More than 2 actins--> more bending forces to calculate
            xap1 = actins[three]->get_xcm() - actins[two]->get_xcm();
            yap1 = actins[three]->get_ycm() - actins[two]->get_ycm();

            Cap1ap1 = xap1*xap1 + yap1*yap1;
            Caap1   = xap1*xa   + yap1*ya;

            theta_ap1 = angle_between_links(two, one); //links[two]->get_angle() - links[one]->get_angle();

            if ( fabs(theta_ap1) < maxSmallAngle )
                coef2 =  kb * ( 1 / sqrt( Caa     * Cap1ap1 ) );
            else
                coef2 =  kb*theta_ap1 / sin(theta_ap1) * ( 1 / sqrt( Caa     * Cap1ap1 ) );

            //Enter loop if more than 3 actins. For 3 actin case, the loop is skipped
            for (unsigned int j = two; j > 1; j--){

                xap2 = actins[j-2]->get_xcm() - actins[j-1]->get_xcm();
                yap2 = actins[j-2]->get_ycm() - actins[j-1]->get_ycm();

                Cap2ap2 = xap2*xap2 + yap2*yap2;
                Cap1ap2 = xap1*xap2 + yap1*yap2;

                theta_ap2 = angle_between_links(j-1, j); //links[j-1]->get_angle() - links[j]->get_angle();
                
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
                
                actins[j]->update_force(forcex, forcey);

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
            
            actins[1]->update_force(forcex, forcey);

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
            
            actins[0]->update_force(fx1, fy1);
        }
    }
}


//wrapper, for fwd_bending_update (and bwd bending update if I ever make it)
void filament::update_bending()
{
    if(links.size() > 1)
        this->fwd_bending_update();

}


int filament::get_nactins(){
    return actins.size();
}

double filament::get_bending_energy(){
    
    double sum = 0, theta;

    if (links.size() < 2) return 0;

    for (unsigned int i = 0; i < links.size() - 1; i++)
    {
        theta = links[i+1]->get_angle() - links[i]->get_angle();
        sum += theta*theta;
    }
    
    return kb*sum/2.0;

}

double filament::get_stretching_energy(){
    
    if(links.size() == 0){
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

double filament::get_kinetic_energy()
{
    return kinetic_energy;
}

double filament::get_potential_energy()
{
    return this->get_stretching_energy() + this->get_bending_energy();
}

double filament::get_total_energy()
{
    return this->get_potential_energy() + this->get_kinetic_energy();
}

void filament::print_thermo()
{
    cout<<"\tKE = "<<this->get_kinetic_energy()<<"\tPE = "<<this->get_potential_energy()<<\
        "\tTE = "<<this->get_total_energy();
}
/*
void filament::lammps_bending_compute()
{
  int i1,i2,i3,n,type;
  double delx1,dely1,delz1,delx2,dely2,delz2;
  double eangle;
  array<double, 3> f1, f3;
  double dtheta,tk;
  double rsq1,rsq2,r1,r2,c,s,a,a11,a12,a22;

  eangle = 0.0;
  
  double **x = atom->x;
  double **f = atom->f;
  int **anglelist = neighbor->anglelist;
  int nanglelist = neighbor->nanglelist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  for (n = 0; n < actins.size()-2; n++) {
    i1 = n;
    i2 = n+1;
    i3 = anglelist[n][2];
    type = anglelist[n][3];

    // 1st bond

    delx1 = x[i1][0] - x[i2][0];
    dely1 = x[i1][1] - x[i2][1];
    delz1 = x[i1][2] - x[i2][2];

    rsq1 = delx1*delx1 + dely1*dely1 + delz1*delz1;
    r1 = sqrt(rsq1);

    // 2nd bond

    delx2 = x[i3][0] - x[i2][0];
    dely2 = x[i3][1] - x[i2][1];
    delz2 = x[i3][2] - x[i2][2];

    rsq2 = delx2*delx2 + dely2*dely2 + delz2*delz2;
    r2 = sqrt(rsq2);

    // angle (cos and sin)

    c = delx1*delx2 + dely1*dely2 + delz1*delz2;
    c /= r1*r2;

    if (c > 1.0) c = 1.0;
    if (c < -1.0) c = -1.0;

    s = sqrt(1.0 - c*c);
    if (s < SMALL) s = SMALL;
    s = 1.0/s;

    // force & energy

    dtheta = acos(c) - theta0[type];
    tk = k[type] * dtheta;

    if (eflag) eangle = tk*dtheta;

    a = -2.0 * tk * s;
    a11 = a*c / rsq1;
    a12 = -a / (r1*r2);
    a22 = a*c / rsq2;

    f1[0] = a11*delx1 + a12*delx2;
    f1[1] = a11*dely1 + a12*dely2;
    f1[2] = a11*delz1 + a12*delz2;
    f3[0] = a22*delx2 + a12*delx1;
    f3[1] = a22*dely2 + a12*dely1;
    f3[2] = a22*delz2 + a12*delz1;

    // apply force to each of 3 atoms

    if (newton_bond || i1 < nlocal) {
      f[i1][0] += f1[0];
      f[i1][1] += f1[1];
      f[i1][2] += f1[2];
    }

    if (newton_bond || i2 < nlocal) {
      f[i2][0] -= f1[0] + f3[0];
      f[i2][1] -= f1[1] + f3[1];
      f[i2][2] -= f1[2] + f3[2];
    }

    if (newton_bond || i3 < nlocal) {
      f[i3][0] += f3[0];
      f[i3][1] += f3[1];
      f[i3][2] += f3[2];
    }

    if (evflag) ev_tally(i1,i2,i3,nlocal,newton_bond,eangle,f1,f3,
                         delx1,dely1,delz1,delx2,dely2,delz2);
  }
}
*/
