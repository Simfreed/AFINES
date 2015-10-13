/*
 * motor.cpp
 *  
 *
 *  Created by Shiladitya Banerjee on 9/3/13.
 *  Copyright 2013 University of Chicago. All rights reserved.
 *
 */

#include "globals.h"
#include "motor.h"
#include "motor_ensemble.h"
#include "filament_ensemble.h"

//motor_ensemble class
template <class filament_ensemble_type>
motor_ensemble<filament_ensemble_type>::motor_ensemble(double mdensity, array<double, 2> myfov, double delta_t, double temp, 
        double mlen, filament_ensemble_type * network, double v0, double stiffness, double ron, double roff, double rend, 
        double actin_len, double vis, vector<array<double,3> > positions, string BC) {
    
    fov = myfov;
    mld =mlen;
    gamma = 0;
    tMove=0;//10;
    f_network=network;
    
    int nm = int(ceil(mdensity*fov[0]*fov[1]));
    cout<<"\nDEBUG: Number of motors:"<<nm<<"\n";

    double alpha = 1, motorx, motory, mang;

    for (int i=0; i< nm; i++) {
        
        if ((unsigned int)i < positions.size()){
            motorx = positions[i][0];
            motory = positions[i][1];
            mang   = positions[i][2];
        }else{
            motorx = rng(-0.5*(fov[0]*alpha-mld),0.5*(fov[0]*alpha-mld));
            motory = rng(-0.5*(fov[1]*alpha-mld),0.5*(fov[1]*alpha-mld));
            mang   = rng(0,2*pi);
        }
        n_motors.push_back(new motor<filament_ensemble_type>( {motorx, motory, mang}, mld, f_network,{0, 0}, {-1,-1}, {-1,-1}, fov, delta_t, temp, 
                    v0, stiffness, ron, roff, rend, actin_len, vis, BC));
        
    }
}

template <class filament_ensemble_type>
motor_ensemble<filament_ensemble_type>::motor_ensemble(vector<vector<double> > motors, array<double, 2> myfov, double delta_t, double temp, 
        double mlen, filament_ensemble_type * network, double v0, double stiffness, double ron, double roff, double rend, 
        double actin_len, double vis, string BC) {
    
    fov = myfov;
    mld =mlen;
    gamma = 0;
    tMove = 0;
    f_network=network;
    
    int nm = motors.size();
    cout<<"\nDEBUG: Number of motors:"<<nm<<"\n";

    double motorx, motory, mang;
    array<int, 2> f_index, l_index, state;

    for (int i=0; i< nm; i++) {
        
        motorx = motors[i][0] + 0.5*motors[i][2];
        motory = motors[i][1] + 0.5*motors[i][3];
        mang   = atan2(motors[i][3], motors[i][2]);
        
        f_index = {int(motors[i][4]), int(motors[i][5])};
        l_index = {int(motors[i][6]), int(motors[i][7])};

        state = {f_index[0] == -1 && l_index[0] == -1 ? 0 : 1, f_index[1] == -1 && l_index[1] == -1 ? 0 : 1};  

        n_motors.push_back(new motor<filament_ensemble_type>( {motorx, motory, mang}, mld, f_network, state, f_index, l_index, fov, delta_t, temp, 
                    v0, stiffness, ron, roff, rend, actin_len, vis, BC));
    }
}

template <class filament_ensemble_type>
motor_ensemble<filament_ensemble_type>::~motor_ensemble( ){ 
    cout<<"DELETING MOTOR ENSEMBLE\n";
    int s = n_motors.size();
    for (int i = 0; i < s; i++){
        delete n_motors[i];
    }
    n_motors.clear();
};

template <class filament_ensemble_type>
int motor_ensemble<filament_ensemble_type>::get_nmotors( ){ 
    return n_motors.size();
}
//check if any motors attached to filaments that no longer exist; 
// if they do, detach them
// Worst case scenario, this is a O(m*n*p) function,
// where m is the number of filaments
//       n is the number of rods per filament
//       p is the number of motors
// However: we don't expect to fracture often, 
// so this loop should rarely if ever be accessed.
    

template <class filament_ensemble_type>
void motor_ensemble<filament_ensemble_type>::check_broken_filaments()
{
    vector<int> broken_filaments = f_network->get_broken();
    array<int, 2> f_index;
    
    for (unsigned int i = 0; i < broken_filaments.size(); i++){
        
        for(unsigned int j = 0; j < n_motors.size(); j++){
            
            f_index = n_motors[j]->get_f_index();

            if(f_index[0] == broken_filaments[i]){
                n_motors[j]->detach_head(0);
                //cout<<"\nDEBUG: detaching head 0 of motor "<<j<<" from broken filament "<<broken_filaments[i];
            }

            if(f_index[1] == broken_filaments[i]){
                n_motors[j]->detach_head(1);
                //cout<<"\nDEBUG: detaching head 1 of motor "<<j<<" from broken filament "<<broken_filaments[i];
            }
        }
    }


}

template <class filament_ensemble_type>
void motor_ensemble<filament_ensemble_type>::motor_walk(double t)
{

    this->check_broken_filaments();
    array<int, 2> s;

    for (unsigned int i=0; i<n_motors.size(); i++) {
       
        s[0]=n_motors[i]->get_states()[0];
        s[1]=n_motors[i]->get_states()[1];

        if (s[0]==0 && s[1]==0) {
            n_motors[i]->attach(t, 0);
            n_motors[i]->attach(t, 1);
        }
        else if (s[0]==0 && s[1]==1) {
            n_motors[i]->attach(t, 0);
        }
        else if (s[0]==1 && s[1]==0) {
            n_motors[i]->attach(t, 1);
        }
    
    }
   

    if (t >= tMove){
        for (unsigned int i=0; i<n_motors.size(); i++) {

            //cout<<"\nDEBUG: motor "<<i;
            s[0]=n_motors[i]->get_states()[0];
            s[1]=n_motors[i]->get_states()[1];
            
            if (s[0]==0 && s[1]==0) {
                n_motors[i]->update_position(t);
            }
            else if (s[0]==0 && s[1]==1) {
                n_motors[i]->update_position(t);
                n_motors[i]->step_onehead(1);
            }
            else if (s[0]==1 && s[1]==0) {
                n_motors[i]->update_position(t);
                n_motors[i]->step_onehead(0);
            }
            else {
                n_motors[i]->update_position(t);
                n_motors[i]->step_twoheads();
            }

            n_motors[i]->actin_update();
        }   
    }
}

template <class filament_ensemble_type>
void motor_ensemble<filament_ensemble_type>::motor_write(ostream& fout)
{
    for (unsigned int i=0; i<n_motors.size(); i++) {
        fout<<"\n"<<n_motors[i]->get_hx()[0]<<"\t"<<n_motors[i]->get_hy()[0]<<"\t"
            <<n_motors[i]->get_hx()[1]-n_motors[i]->get_hx()[0]<<"\t"
            <<n_motors[i]->get_hy()[1]-n_motors[i]->get_hy()[0]<<"\t"
            <<n_motors[i]->get_f_index()[0]<<"\t"<<n_motors[i]->get_f_index()[1]<<"\t"
            <<n_motors[i]->get_l_index()[0]<<"\t"<<n_motors[i]->get_l_index()[1];
    } 
}

template <class filament_ensemble_type>
void motor_ensemble<filament_ensemble_type>::add_motor(motor<filament_ensemble_type> * m)
{
    n_motors.push_back(m);
}

template <class filament_ensemble_type>
void motor_ensemble<filament_ensemble_type>::set_shear(double g)
{
    for (unsigned int i=0; i<n_motors.size(); i++)
        n_motors[i]->set_shear(g);
    
    gamma = g;
}

template class motor_ensemble<ATfilament_ensemble>;
template class motor_ensemble<lammps_filament_ensemble>;
template class motor_ensemble<baoab_filament_ensemble>;
template class motor_ensemble<langevin_leapfrog_filament_ensemble>;
