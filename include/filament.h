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

        filament(array<double, 2> myfov, array<int, 2> mynq, double deltat, double temp, double shear, 
                double frac, double bending_stiffness, string bndcnd);
        
        filament(array<double, 3> startpos, int nactin, array<double,2> myfov, array<int,2> mynq,
                double vis, double deltat, double temp, bool isStraight,
                double actinLength, double linkLength, double stretching, double bending, double fracture, string bc); 

        filament(vector<actin *> actinvec, array<double, 2> myfov, array<int, 2> mynq, double linkLength, double stretching_stiffness, double bending_stiffness, 
                double deltat, double temp, double fracture, double gamma, string bc);
       
        filament();
        
        ~filament();
    
        void set_shear(double g);

        void update_shear(); 
        
        vector<filament *> update_stretching();
        
        void update_bending();
        
        void update_positions(double t);
        
        array<double, 2> boundary_check(int i, double t, double vx, double vy);
        
        void update(double t);
        
        actin * get_actin(int i);
        
        Link * get_link(int i);

        int get_nlinks();

        vector<vector<array<int, 2> > > get_quadrants();
        
        string write_actins();
        
        string write_links();
        
        string to_string();
        
        string write_thermo();
        
        vector<actin *> get_actins(unsigned int first, unsigned int last);
        
        vector<filament *> fracture(int node);
        
        void update_forces(int index, double f1, double f2);
        
        bool operator==(const filament& that);
        
        void add_actin(actin * a, double l0, double kl);
        
        void set_BC(string s);

        string get_BC();
       
        inline double angle_between_links(int i, int j);

        void fwd_bending_update();
        
        void bwd_bending_update();
        
        int get_nactins();

        double get_bending_energy();

        double get_stretching_energy();

        double get_kinetic_energy();
        
        double get_potential_energy();
        
        double get_total_energy();
    
        void print_thermo();

    protected:
        
        double kb, gamma, temperature, dt, fracture_force, kinetic_energy;
        
        array<double,2> fov;
        array<int,2> nq;
        vector<array<double, 2> > prv_rnds;
        vector<actin *> actins;
        vector<Link *> links;
        string BC;
};

class baoab_filament : public filament
{

    //using filament::filament;
    public: 
        
        baoab_filament(array<double, 3> startpos, int nactin, array<double,2> myfov, array<int,2> mynq,
                double vis, double deltat, double temp, bool isStraight,
                double actinLength, double linkLength, double stretching, double bending, double fracture, string bc); 

        baoab_filament(vector<actin *> actinvec, array<double, 2> myfov, array<int, 2> mynq, double linkLength, double stretching_stiffness, double bending_stiffness, 
                double deltat, double temp, double fracture, double gamma, string bc);
        
        ~baoab_filament();

        void update_velocities_O(double t);

        void update_velocities_B();

        void update_positions(double t);
        
        vector<baoab_filament *> update_stretching();
        
        vector<baoab_filament *> fracture(int node);
        
        bool operator==(const baoab_filament& that);

    protected:
        double a, b, mass;
};

class lammps_filament : public filament
{

    //using filament::filament;
    public:
       
        lammps_filament(array<double, 3> startpos, int nactin, array<double,2> myfov, array<int,2> mynq,
                double vis, double deltat, double temp, bool isStraight,
                double actinLength, double linkLength, double stretching, double bending, double fracture, string bc); 

        lammps_filament(vector<actin *> actinvec, array<double, 2> myfov, array<int, 2> mynq, double linkLength, double stretching_stiffness, double bending_stiffness, 
                double deltat, double temp, double fracture, double gamma, string bc);
        
        lammps_filament(array<double, 2> myfov, array<int, 2> mynq, double deltat, double temp, double shear, 
                double frac, double bending_stiffness, string bndcnd);
       
        lammps_filament();
        
        ~lammps_filament();
        
        void set_mass(double m);

        void update_brownian();

        void update_drag();

        void update_positions(double t);
        
        vector<lammps_filament *> update_stretching();
    
        vector<lammps_filament *> fracture(int node);

        bool operator==(const lammps_filament& that);
    
    protected:
        double mass;
};       

class langevin_leapfrog_filament : public filament
{

    //using filament::filament;
    public:
       
        langevin_leapfrog_filament(array<double, 3> startpos, int nactin, array<double,2> myfov, array<int,2> mynq,
                double vis, double deltat, double temp, bool isStraight,
                double actinLength, double linkLength, double stretching, double bending, double fracture, string bc); 

        langevin_leapfrog_filament(vector<actin *> actinvec, array<double, 2> myfov, array<int, 2> mynq, double linkLength, double stretching_stiffness, double bending_stiffness, 
                double deltat, double temp, double fracture, double gamma, string bc);
        
        langevin_leapfrog_filament(array<double, 2> myfov, array<int, 2> mynq, double deltat, double temp, double shear, 
                double frac, double bending_stiffness, string bndcnd);
       
        langevin_leapfrog_filament();
        
        ~langevin_leapfrog_filament();
        
        void set_mass(double m);

        void reset_velocity();

        void update_velocity_brownian();

        void update_velocity_drag();

        void update_velocity_int_forces();
        
        void update_positions(double t);
        
        vector<langevin_leapfrog_filament *> update_stretching();
    
        vector<langevin_leapfrog_filament *> fracture(int node);

        bool operator==(const langevin_leapfrog_filament& that);
    
    protected:
        double mass, a, b, c;
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

