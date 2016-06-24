/*
 *  filament.cpp
 *  
 *
 *  Created by Simon Freedman and Shiladitya Banerjee
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
    fracture_force = 1000000;
    BC = "REFLECTIVE";
    kinetic_energy = 0;
    gamma = 0;
    delrx=0;
    damp = infty;
    y_thresh=2;
    kToverLp=0;
    bd_prefactor = sqrt(temperature/(2*dt*damp));
}

filament::filament(array<double, 2> myfov, array<int, 2> mynq, double deltat, double temp, double shear, 
        double frac, double bending_stiffness, string bndcnd)
{
    fov             = myfov;
    nq              = mynq;
    dt              = deltat;
    temperature     = temp;
    gamma           = shear;
    delrx           = 0;
    fracture_force  = frac;
    kb              = bending_stiffness;
    BC              = bndcnd;
    kinetic_energy  = 0;
    damp            = infty;
    y_thresh        = 1;
    kToverLp        = temperature * temperature / (bending_stiffness*1); // 1 is the link length, here
    bd_prefactor = sqrt(temperature/(2*dt*damp));

}

filament::filament(array<double, 3> startpos, int nactin, array<double, 2> myfov, array<int, 2> mynq, double visc, 
        double deltat, double temp, bool isStraight, double actinRadius, double linkLength, double stretching_stiffness,
        double max_ext_ratio, double bending_stiffness, double frac_force, string bdcnd)
{
    
    fov = myfov;
    nq = mynq;
    dt = deltat;
    temperature = temp;
    gamma = 0;
    delrx = 0;
    fracture_force = frac_force;
    BC = bdcnd;
    kb = bending_stiffness;
    kinetic_energy  = 0;
    
    damp = 4*pi*actinRadius*visc;
    y_thresh = 1;
    
    kToverLp = temperature * temperature / (kb * linkLength); 
    bd_prefactor = sqrt(temperature/(2*dt*damp));

    double phi, variance;
    array<double, 2> next_pos;
    //the start of the polymer: 
    actins.push_back(new actin( startpos[0], startpos[1], actinRadius, visc));
    prv_rnds.push_back({0,0});
    phi = startpos[2];
    
    if (temp != 0) variance = temp/bending_stiffness;
    else variance = 0;

    for (int j = 1; j < nactin; j++) {

        //next_pos = boundary_check(j-1, actins[j-1]->get_xcm() + linkLength*cos(phi), actins[j-1]->get_ycm() + linkLength*sin(phi));
        next_pos = pos_bc(BC, delrx, dt, fov, 
                {linkLength*cos(phi)/dt, linkLength*sin(phi)/dt},
                {actins[j-1]->get_xcm() + linkLength*cos(phi), actins[j-1]->get_ycm() + linkLength*sin(phi)});
        actins.push_back( new actin(next_pos[0], next_pos[1], actinRadius, visc) );
        prv_rnds.push_back({0,0});
        links.push_back( new Link(linkLength, stretching_stiffness, max_ext_ratio, this, {j-1, j}, fov, nq) );  
        links[j-1]->step(BC, delrx);  
        
        // Calculate the Next angle on the actin polymer
        if (!isStraight) phi += rng_n(0, variance);
    
    }
   
}

filament::filament(vector<actin *> actinvec, array<double, 2> myfov, array<int, 2> mynq, double linkLength, 
        double stretching_stiffness, double max_ext_ratio, double bending_stiffness, 
        double deltat, double temp, double frac_force, double g, string bdcnd)
{
    
    dt = deltat;
    temperature = temp;
    fracture_force = frac_force;
    gamma = g; 
    delrx = 0;
    BC = bdcnd;
    kb = bending_stiffness;
    fov = myfov;
    nq = mynq;
    y_thresh = 1;

    kToverLp = temperature * temperature / (kb * linkLength); 

    if (actinvec.size() > 0)
    {
        actins.push_back(new actin(*(actinvec[0])));
        prv_rnds.push_back({0,0});
        damp = actins[0]->get_friction();
    }
    
    //Link em up
    if (actinvec.size() > 1){
        for (unsigned int j = 1; j < actinvec.size(); j++) {

            actins.push_back(new actin(*(actinvec[j])));
            links.push_back( new Link(linkLength, stretching_stiffness, max_ext_ratio, this, {(int)j-1, (int)j}, fov, nq) );  
            links[j-1]->step(BC, delrx);
            prv_rnds.push_back({0,0});
            
        }
    }
    
    bd_prefactor = sqrt(temperature/(2*dt*damp));

}

filament::~filament(){
    
    //cout<<"DELETING FILAMENT\n";
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
    prv_rnds.clear();
}

void filament::add_actin(actin * a, double linkLength, double stretching_stiffness, double max_ext_ratio){
    
    actins.push_back(new actin(*a));
    prv_rnds.push_back({0,0});    
    if (actins.size() > 1){
        int j = (int) actins.size() - 1;
        links.push_back( new Link(linkLength, stretching_stiffness, max_ext_ratio, this, {j-1,  j}, fov, nq ) );  
        links[j-1]->step(BC, delrx);
    }
    if (damp == infty)
        damp = a->get_friction();
}

vector<vector<array<int,2> > > filament::get_quadrants()
{
    //should return a map between actin and x, y coords of quadrant
    vector<vector<array<int,2> > > quads;
    for (unsigned int i=0; i < links.size(); i++){ 
        links[i]->quad_update(BC, delrx);
        quads.push_back(links[i]->get_quadrants());
    }
    
    return quads;
}
/*
multimap<int, array<int, 2> > filament::get_quadrants()
{
    //should return a map between actin and x, y coords of quadrant
    multipmap<int, array<int, 2> > quads;
    for (unsigned int i=0; i < links.size(); i++){ 
        links[i]->quad_update(BC, delrx);
        quads.push_back(links[i]->get_quadrants());
    }
    
    return quads;
}*/

void filament::set_y_thresh(double y){
    y_thresh = y;
}

void filament::update_positions()
{
    double vx, vy;
    array<double, 2> new_rnds;
    array<double, 2> newpos;
    kinetic_energy = 0;  
    double top_y = y_thresh*fov[1]/2.; 
    int sa = int(actins.size());
    int la = int(links.size());
    for (int i = 0; i < sa; i++){

        if (fabs(actins[i]->get_ycm()) > top_y) continue;
     
        new_rnds = {rng_n(), rng_n()};
        vx  = (actins[i]->get_force()[0])/damp  + bd_prefactor*(new_rnds[0] + prv_rnds[i][0]);
        vy  = (actins[i]->get_force()[1])/damp  + bd_prefactor*(new_rnds[1] + prv_rnds[i][1]);
//        cout<<"\nDEBUG: Fx("<<i<<") = "<<actins[i]->get_force()[0]<<"; v = ("<<vx<<" , "<<vy<<")";
       
        prv_rnds[i] = new_rnds;
        //cout<<"\nDEBUG: actin force = ("<<actins[i]->get_force()[0]<<" , "<<actins[i]->get_force()[1]<<")";
        kinetic_energy += vx*vx + vy*vy;
        newpos = pos_bc(BC, delrx, dt, fov, {vx, vy}, {actins[i]->get_xcm() + vx*dt, actins[i]->get_ycm() + vy*dt});
        actins[i]->set_xcm(newpos[0]);
        actins[i]->set_ycm(newpos[1]);
        actins[i]->reset_force(); 
    }

    for (int i = 0; i < la; i++)
        links[i]->step(BC, delrx);

}

void filament::update_positions_range(int lo, int hi)
{
    double vx, vy;
    array<double, 2> new_rnds;
    array<double, 2> newpos;
    kinetic_energy = 0;  
    double top_y = y_thresh*fov[1]/2.; 

    int low = max(0, lo);
    int high = min(hi, (int)actins.size());

    for (int i = low; i < high; i++){
       
        if (fabs(actins[i]->get_ycm()) > top_y) continue;
     
        new_rnds = {rng_n(), rng_n()};
        vx  = (actins[i]->get_force()[0])/damp  + bd_prefactor*(new_rnds[0] + prv_rnds[i][0]);
        vy  = (actins[i]->get_force()[1])/damp  + bd_prefactor*(new_rnds[1] + prv_rnds[i][1]);
//        cout<<"\nDEBUG: Fx("<<i<<") = "<<actins[i]->get_force()[0]<<"; v = ("<<vx<<" , "<<vy<<")";
       
        prv_rnds[i] = new_rnds;
        //cout<<"\nDEBUG: actin force = ("<<actins[i]->get_force()[0]<<" , "<<actins[i]->get_force()[1]<<")";
        kinetic_energy += vx*vx + vy*vy;
        //newpos = boundary_check(i, actins[i]->get_xcm() + vx*dt, actins[i]->get_ycm() + vy*dt); 
        newpos = pos_bc(BC, delrx, dt, fov, {vx, vy}, {actins[i]->get_xcm() + vx*dt, actins[i]->get_ycm() + vy*dt});
        actins[i]->set_xcm(newpos[0]);
        actins[i]->set_ycm(newpos[1]);
        actins[i]->reset_force(); 
    }

    for (unsigned int i = 0; i < links.size(); i++)
        links[i]->step(BC, delrx);

}

vector<filament *> filament::update_stretching(double t)
{
    vector<filament *> newfilaments;

    if(links.size() == 0)
        return newfilaments;
   
    for (unsigned int i=0; i < links.size(); i++) {
        links[i]->update_force(BC, delrx);
        //links[i]->update_force_fraenkel_fene(BC, delrx);
        if (hypot(links[i]->get_force()[0], links[i]->get_force()[1]) > fracture_force){
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

void filament::update_delrx(double shear_dist){
    delrx = shear_dist;
}

void filament::update_shear(double t){
    
    double local_shear;
    for (unsigned int i = 0; i < actins.size(); i++){
        local_shear = delrx * actins[i]->get_ycm() / fov[1];
        actins[i]->set_xcm(actins[i]->get_xcm() + local_shear);
        //cout<<"\nDEBUG: local_shear = "<<local_shear;
    }
}

void filament::update_d_strain(double g){
    
    for (unsigned int i = 0; i < actins.size(); i++){
        actins[i]->set_xcm(actins[i]->get_xcm() + g * actins[i]->get_ycm() / fov[1]);
    }
}

void filament::update_forces(int index, double f1, double f2)
{
    actins[index]->update_force(f1,f2);
}

void filament::pull_on_ends(double f)
{
    if (actins.size() < 2) return;
    int last = actins.size() - 1; 
    array<double, 2> dr = rij_bc(BC, actins[last]->get_xcm() - actins[0]->get_xcm(),  
                                     actins[last]->get_ycm() - actins[0]->get_ycm(), fov[0], fov[1], delrx);
    double ang = atan2( dr[1], dr[0]);
    
    actins[ 0  ]->update_force(-0.5*f*cos(ang), -0.5*f*sin(ang));
    actins[last]->update_force( 0.5*f*cos(ang),  0.5*f*sin(ang));
}

void filament::affine_pull(double f)
{
    if (actins.size() < 2) return;
    int last = actins.size() - 1; 
    array<double, 2> dr = rij_bc(BC, actins[last]->get_xcm() - actins[0]->get_xcm(),  
                                     actins[last]->get_ycm() - actins[0]->get_ycm(), fov[0], fov[1], delrx);
    double ang = atan2( dr[1], dr[0]);
    //cout<<"\nDEBUG: angle = "<<ang;
    double frac, fcos = f*cos(ang), fsin = f*sin(ang);

    for (int i = 0; i <= last; i++){
        frac = (double(i)/double(last)-0.5);
        actins[i]->update_force(frac*fcos, frac*fsin);
    }
}

void filament::set_shear(double g){
    gamma = g;
    max_shear = gamma*fov[1]*0.5;
}

string filament::write_actins(int fil){
    string all_actins;
    for (unsigned int i =0; i < actins.size(); i++)
    {
        all_actins += actins[i]->write() + "\t" + std::to_string(fil);
    }

    return all_actins;
}

string filament::write_links(int fil){
    string all_links;
    for (unsigned int i =0; i < links.size(); i++)
    {
        all_links += links[i]->write(BC, delrx) + "\t" + std::to_string(fil);
    }

    return all_links;
}

string filament::write_thermo(int fil)
{
    return "\n" + std::to_string(this->get_kinetic_energy()) + \
        "\t" + std::to_string(this->get_potential_energy()) + \
        "\t" + std::to_string(this->get_total_energy()) + "\t" + std::to_string(fil);
}

vector<actin *> filament::get_actins(unsigned int first, unsigned int last)
{
    vector<actin *> newactins;
    for (unsigned int i = first; i < last; i++)
    {
        if (i >= actins.size())
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
                new filament(lower_half, fov, nq, links[0]->get_length(), links[0]->get_kl(), links[0]->get_fene_ext(), kb, 
                    dt, temperature, fracture_force, gamma, BC));
    if (upper_half.size() > 0)
        newfilaments.push_back(
                new filament(upper_half, fov, nq, links[0]->get_length(), links[0]->get_kl(), links[0]->get_fene_ext(), kb, 
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
    double theta = links[i]->get_angle() - links[j]->get_angle();
    return theta - 2*pi*floor(theta/(2*pi)+0.5); // Keep angles between -Pi and Pi
}

void filament::fwd_bending_update(double t)
{
    //Calculate the force at each Link position as outlined by 
    // Allen and Tildesley in Computer Simulation of Liquids, Appendix C
    // Here we let db = (xb, yb) = (linkx[b] - linkx[b-1], linky[b] - linky[b-1])
 
    array<double, 2> ram1, ra, rap1, rap2;
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
       
        ram1 = links[0]->get_disp();
        ra = links[1]->get_disp();

        Cam1am1 = ram1[0]*ram1[0] + ram1[1]*ram1[1];
        Caa     = ra[0]*ra[0]   + ra[1]*ra[1];
        Caam1   = ram1[0]*ra[0]   + ram1[1]*ra[1];

        theta_a   = angle_between_links(1,0);

        if( fabs(theta_a) < maxSmallAngle )
            coef1 = -kb * ( 1 / sqrt( Cam1am1 * Caa     ) ); 
        else
            coef1 = -kb*theta_a   / sin(theta_a)   * ( 1 / sqrt( Cam1am1 * Caa     ) ); 

        //2 actins --> bending force on last link; only has one term
        if (actins.size() == 3)
        {

            fx1 = coef1 * ( Caam1/Caa * ra[0] - ram1[0] );
            fy1 = coef1 * ( Caam1/Caa * ra[1] - ram1[1] );
            actins[2]->update_force(fx1, fy1);

        }
        else
        {
            //More than 2 actins--> more bending forces to calculate
            rap1 = links[2]->get_disp();

            Cap1ap1 = rap1[0]*rap1[0] + rap1[1]*rap1[1];
            Caap1   = rap1[0]*ra[0]   + rap1[1]*ra[1];

            theta_ap1 = angle_between_links(2,1); //links[2]->get_angle() - links[1]->get_angle();
            
            if ( fabs(theta_ap1) < maxSmallAngle )
                coef2 =  kb * ( 1 / sqrt( Caa     * Cap1ap1 ) );
            else
                coef2 =  kb*theta_ap1 / sin(theta_ap1) * ( 1 / sqrt( Caa     * Cap1ap1 ) );

            //Enter loop if more than 3 actins. For 3 actin case, the loop is skipped
            for (unsigned int j = 2; j < actins.size() - 2; j++){

                
                rap2 = links[j+1]->get_disp();
                
                Cap2ap2 = rap2[0]*rap2[0] + rap2[1]*rap2[1];
                Cap1ap2 = rap1[0]*rap2[0] + rap1[1]*rap2[1];

                theta_ap2 = angle_between_links(j+1,j); //links[j+1]->get_angle() - links[j]->get_angle();
                
                if ( fabs(theta_ap2) < maxSmallAngle )
                    coef3 = -kb * ( 1 / sqrt( Cap1ap1 * Cap2ap2 ) );
                else
                    coef3 = -kb*theta_ap2 / sin(theta_ap2) * ( 1 / sqrt( Cap1ap1 * Cap2ap2 ) );
                
                fx1 = coef1 * ( Caam1/Caa * ra[0] - ram1[0] );
                fy1 = coef1 * ( Caam1/Caa * ra[1] - ram1[1] );

                fx2 = coef2 * ( (1 + Caap1/Cap1ap1) * rap1[0] - (1 + Caap1/Caa) * ra[0] ); 
                fy2 = coef2 * ( (1 + Caap1/Cap1ap1) * rap1[1] - (1 + Caap1/Caa) * ra[1] ); 

                fx3 = coef3 * ( rap2[0] - Cap1ap2/Cap1ap1 * rap1[0]);
                fy3 = coef3 * ( rap2[1] - Cap1ap2/Cap1ap1 * rap1[1]);

                forcex = fx1 + fx2 + fx3;
                forcey = fy1 + fy2 + fy3;
                
                actins[j]->update_force(forcex, forcey);

                //increment all variables for next iteration:
                ram1 = ra;
                ra   = rap1;
                rap1 = rap2;

                Cam1am1 = Caa;
                Caa     = Cap1ap1;
                Cap1ap1 = Cap2ap2;

                Caam1 = Caap1;
                Caap1 = Cap1ap2;

                coef1 = -1*coef2;
                coef2 = -1*coef3;

            }

            //LAST TWO actinS ON THE FILAMENT:
            fx1 = coef1 * ( Caam1/Caa * ra[0] - ram1[0] );
            fy1 = coef1 * ( Caam1/Caa * ra[1] - ram1[1] );
            fx2 = coef2 * ( (1 + Caap1/Cap1ap1) * rap1[0] - (1 + Caap1/Caa) * ra[0] ); 
            fy2 = coef2 * ( (1 + Caap1/Cap1ap1) * rap1[1] - (1 + Caap1/Caa) * ra[1] ); 
            
            //cout<<"\nDEBUG: 2nd to last actin: (fx1, fy1) + (fx2, fy2) = ( "<<fx1<< " , "<<fy1<<" ) + ( "<<fx2<<" , " <<fy2<<" )";
            
            forcex = fx1 + fx2;
            forcey = fy1 + fy2;
            
            actins[ actins.size() - 2 ]->update_force(forcex, forcey);

            /*INCREMENT*/
            ram1 = ra;
            ra   = rap1;

            Caa   = Cap1ap1;
            Caam1 = Caap1;

            coef1 = -1*coef2;

            fx1 = coef1 * ( Caam1/Caa * ra[0] - ram1[0] );
            fy1 = coef1 * ( Caam1/Caa * ra[1] - ram1[1] );
            
            actins[ actins.size() - 1 ]->update_force(fx1, fy1);
        }
    }
}

void filament::bwd_bending_update(double t)
{
    
    //same thing as fwd_bending_update, but starting from the other end of the filament
    array<double, 2> ram1, ra, rap1, rap2; 
    double Cam1am1, Caa, Cap1ap1, Cap2ap2, Caam1, Caap1, Cap1ap2;
    double theta_a, theta_ap1, theta_ap2;
    double coef1, coef2, coef3;
    double fx1, fx2, fx3, fy1, fy2, fy3;
    
    double forcex, forcey;
    //initialize all NODE forces to be 0
   
    int one = links.size() - 1, two = links.size()-2, three = links.size() - 3;
    //single actin --> no bending energy 
    if (actins.size() > 2)
    { 
        //First two actins won't have any bending forces: 
        
        // More than 1 actin--> bending forces to calculate
        ram1 = links[one]->get_neg_disp(); 
        ra = links[two]->get_neg_disp(); 

        Cam1am1 = ram1[0]*ram1[0] + ram1[1]*ram1[1];
        Caa     = ra[0]  *ra[0]   + ra[1]  *ra[1];
        Caam1   = ram1[0]*ra[0]   + ram1[1]*ra[1];

        theta_a   = angle_between_links(two, one); //links[one]->get_angle() - links[zero]->get_angle();

        if( fabs(theta_a) < maxSmallAngle )
            coef1 = -kb * ( 1 / sqrt( Cam1am1 * Caa     ) ); 
        else
            coef1 = -kb*theta_a   / sin(theta_a)   * ( 1 / sqrt( Cam1am1 * Caa     ) ); 

        //2 actins --> bending force on last link; only has one term
        if (actins.size() == 3)
        {

            fx1 = coef1 * ( Caam1/Caa * ra[0] - ram1[0] );
            fy1 = coef1 * ( Caam1/Caa * ra[1] - ram1[1] );
            //cout<<"\nDEBUG: magnitude of bending forces at actins: ( "<<fx1<<" , "<<fy1<<" )";
            
            actins[two]->update_force(fx1, fy1);

        }
        else
        {
            //More than 2 actins--> more bending forces to calculate
            rap1 = links[three]->get_neg_disp();
            Cap1ap1 = rap1[0]*rap1[0] + rap1[1]*rap1[1];
            Caap1   = rap1[0]*ra[0]   + rap1[1]*ra[1];

            theta_ap1 = angle_between_links(three, two); //links[two]->get_angle() - links[one]->get_angle();

            if ( fabs(theta_ap1) < maxSmallAngle )
                coef2 =  kb * ( 1 / sqrt( Caa     * Cap1ap1 ) );
            else
                coef2 =  kb*theta_ap1 / sin(theta_ap1) * ( 1 / sqrt( Caa     * Cap1ap1 ) );

            //Enter loop if more than 3 actins. For 3 actin case, the loop is skipped
            for (unsigned int j = two; j > 1; j--){

                rap2 = links[j-2]->get_neg_disp();
                Cap2ap2 = rap2[0]*rap2[0] + rap2[1]*rap2[1];
                Cap1ap2 = rap1[0]*rap2[0] + rap1[1]*rap2[1];

                theta_ap2 = angle_between_links(j-2, j-1); //links[j-1]->get_angle() - links[j]->get_angle();
                
                if ( fabs(theta_ap2) < maxSmallAngle )
                    coef3 = -kb * ( 1 / sqrt( Cap1ap1 * Cap2ap2 ) );
                else
                    coef3 = -kb*theta_ap2 / sin(theta_ap2) * ( 1 / sqrt( Cap1ap1 * Cap2ap2 ) );
                
                fx1 = coef1 * ( Caam1/Caa * ra[0] - ram1[0] );
                fy1 = coef1 * ( Caam1/Caa * ra[1] - ram1[1] );

                fx2 = coef2 * ( (1 + Caap1/Cap1ap1) * rap1[0] - (1 + Caap1/Caa) * ra[0] ); 
                fy2 = coef2 * ( (1 + Caap1/Cap1ap1) * rap1[1] - (1 + Caap1/Caa) * ra[1] ); 

                fx3 = coef3 * ( rap2[0] - Cap1ap2/Cap1ap1 * rap1[0]);
                fy3 = coef3 * ( rap2[1] - Cap1ap2/Cap1ap1 * rap1[1]);

                forcex = fx1 + fx2 + fx3;
                forcey = fy1 + fy2 + fy3;
                
                actins[j]->update_force(forcex, forcey);

                //increment all variables for next iteration:
                ram1 = ra;
                ra   = rap1;
                rap1 = rap2;

                Cam1am1 = Caa;
                Caa     = Cap1ap1;
                Cap1ap1 = Cap2ap2;

                Caam1 = Caap1;
                Caap1 = Cap1ap2;

                coef1 = -1*coef2;
                coef2 = -1*coef3;

            }

            //LAST TWO actinS ON THE FILAMENT:
            fx1 = coef1 * ( Caam1/Caa * ra[0] - ram1[0] );
            fy1 = coef1 * ( Caam1/Caa * ra[1] - ram1[1] );
            fx2 = coef2 * ( (1 + Caap1/Cap1ap1) * rap1[0] - (1 + Caap1/Caa) * ra[0] ); 
            fy2 = coef2 * ( (1 + Caap1/Cap1ap1) * rap1[1] - (1 + Caap1/Caa) * ra[1] ); 
            
            //cout<<"\nDEBUG: 2nd to last actin: (fx1, fy1) + (fx2, fy2) = ( "<<fx1<< " , "<<fy1<<" ) + ( "<<fx2<<" , " <<fy2<<" )";
            
            forcex = fx1 + fx2;
            forcey = fy1 + fy2;
            actins[1]->update_force(forcex, forcey);

            /*INCREMENT*/
            ram1 = ra;
            ra   = rap1;

            Caa   = Cap1ap1;
            Caam1 = Caap1;

            coef1 = -1*coef2;

            fx1 = coef1 * ( Caam1/Caa * ra[0] - ram1[0] );
            fy1 = coef1 * ( Caam1/Caa * ra[1] - ram1[1] );
            //cout<<"\nDEBUG: Last actin: (fx1, fy1) = ( "<<fx1<< " , "<<fy1<<" )";

            actins[0]->update_force(fx1, fy1);
        }
    }
}


//wrapper, for fwd_bending_update (and bwd bending update if I ever make it)
void filament::update_bending(double t)
{
    if(links.size() > 1){
        this->fwd_bending_update(t);
        this->bwd_bending_update(t);
    }
}

int filament::get_nactins(){
    return actins.size();
}

int filament::get_nlinks(){
    return links.size();
}

double filament::get_bending_energy(){
    
    double sum = 0, theta;

    if (links.size() < 2) return 0;

    for (unsigned int i = 0; i < links.size() - 1; i++)
    {
        theta = angle_between_links(i+1, i);//links[i+1]->get_angle() - links[i]->get_angle();
        sum += theta*theta;
    }
    
    return kb*sum; //doubled because of backward and forward bending updates

}

double filament::get_stretching_energy()
{
    
    double u = 0;
    for (unsigned int i = 0; i < links.size(); i++)
        //u += links[i]->get_stretching_energy_fene(BC, delrx);
        u += links[i]->get_stretching_energy();
    
    return u;

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

array<double, 2> filament::get_bead_position(int n)
{
    return {actins[n]->get_xcm(), actins[n]->get_ycm()};
}

void filament::print_thermo()
{
    cout<<"\tKE = "<<this->get_kinetic_energy()<<"\tPE = "<<this->get_potential_energy()<<\
        "\tTE = "<<this->get_total_energy();
}

double filament::get_end2end()
{
    if (actins.size() < 2) 
        return 0;
    else 
        return dist_bc(BC, actins[actins.size() - 1]->get_xcm() - actins[0]->get_xcm(),  
                           actins[actins.size() - 1]->get_ycm() - actins[0]->get_ycm(), fov[0], fov[1], delrx);
} 
