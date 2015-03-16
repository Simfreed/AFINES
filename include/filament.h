/*
 *  filament.cpp
 *  
 *
 *  Created by Simon Freedman on 12/22/2014
 *  Copyright 2014 University of Chicago. All rights reserved.
 *
 */

//=====================================
//include guard
#ifndef __FILAMENT_H_INCLUDED__
#define __FILAMENT_H_INCLUDED__

//=====================================
// forward declared dependencies

//=====================================
//included dependences
#include "string"
#include "vector"
#include "actin.h"
#include "Link.h"
#include "globals.h"

//=====================================
//actin filament class
class filament
{
    public:

        filament(array<double, 3> startpos, int nactin, array<double,2> myfov, array<int,2> mynq,
                double vis, double deltat, double temp, bool isStraight,
                double actinLength, double linkLength, double stretching, double bending, double fracture, string bc); 

        filament(vector<actin *> actinvec, array<double, 2> myfov, array<int, 2> mynq, double linkLength, double stretching_stiffness, double bending_stiffness, 
                double deltat, double temp, double fracture, double gamma, string bc);
       
        filament();
        
        ~filament();
    
        void set_shear(double g);

        void update(double t);
        
        void update_bending();
        
        vector<filament *> update_stretching();

        void update_shear(); 
        
        actin * get_actin(int i);
        
        Link * get_link(int i);

        int get_nlinks();

        vector<vector<array<int, 2> > > get_quadrants();
        
        string write_actins();
        
        string write_links();
        
        string to_string();
        
        vector<actin *> get_actins(unsigned int first, unsigned int last);
        
        vector<filament *> fracture(int node);
        
        void update_forces(int index, double f1, double f2);
        
        bool operator==(const filament& that);
        
        void add_actin(actin * a, double l0, double kl);
        
        void set_BC(string s);

        string get_BC();
       
        void fwd_bending_update();
        
        void bwd_bending_update();
        
        int get_nactins();

        double get_bending_energy();

        double get_stretching_energy();

        double get_total_energy();
    
    protected:
        
        double kb, gamma, temperature, dt, fracture_force;
        
        array<double,2> fov;
        array<int,2> nq;
        vector<actin *> actins;
        vector<Link *> links;
        string BC;
};

// Filament class that is closer to the Nedelec and Foethke model than the above one

class NFfilament : public filament
{
    public:
        
        vector<double> get_P_matrix();

        vector<double> get_A_matrix();

        vector<double> get_B_matrix();

        vector<double> get_G_matrix();

        double* get_mobility();
        
        vector<double *> get_mobility_matrix();

        double get_tau();
        
        vector<NFfilament *> update_stretching();

        void update(double t);
        
        map<array<int,2>, double> get_J_matrix();

};

#endif

