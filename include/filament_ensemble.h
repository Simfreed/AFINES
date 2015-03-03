/*
 *  filament_ensemble.h
 *  
 *
 *  Authors : Shiladitya Banerjee, Simon Freedman
 *  Copyright 2013 University of Chicago. All rights reserved.
 *
 */

//=====================================
//include guard
#ifndef __FILAMENT_ENSEMBLE_H_INCLUDED__
#define __FILAMENT_ENSEMBLE_H_INCLUDED__

//=====================================
// forward declared dependencies

//=====================================
//included dependences
#include "string"
#include "vector"
#include "map"
#include "filament.h"

//=====================================
//filament network class
template <class filament_type>
class filament_ensemble
{
    public:
        filament_ensemble();

        ~filament_ensemble();
        
        void quad_update_monomer(int i);
        
        void quad_update();

        vector<filament_type *> * get_network();

        map<vector<int>, double> get_dist(double x, double y);

        array<double,2> get_direction(int fil, int rod);

        array<double,2> get_intpoints(int fil, int rod, double xp, double yp);

        array<double,2> get_start(int fil, int rod);
        
        array<double,2> get_end(int fil, int rod);
        
        array<double,3> get_forces(int fil, int rod);
        
        double get_int_direction(int fil, int rod, double xp, double yp);

        double get_xcm(int fil, int rod);
        
        double get_ycm(int fil, int rod);

        double get_angle(int fil, int rod);

        double get_alength(int fil, int rod);
        
        void update(double t);

        void update_forces(int fil, int rod, double f1, double f2, double f3);

        void write_rods(ofstream& fout);
        
        void write_links(ofstream& fout);
        
        void set_straight_filaments(bool is_straight);

        void set_shear_rate(double);

        void update_shear();
        
        bool is_polymer_start(int f, int r);

        void set_fov(double x, double y);

        void set_nq(double x, double y);

        void set_ld(double l);

        void set_visc(double v);

        vector<int> get_broken();

        void clear_broken();
        
    protected:

        double dt, temperature, rho, ld, link_ld, visc, gamma;
        int npolymer, nmon;
        bool straight_filaments = false;
        
        array<double,2> fov, view;
        array<int, 2> nq;
        vector<int> broken_filaments, empty_vector;
        
        map<int, map<int,vector< vector<int> > > > quad_fils;
        map<vector<int>,double> t_map;
    
        vector<filament_type *> network;
};

class ATfilament_ensemble:
    public filament_ensemble<filament>
{

    public:
        
        ATfilament_ensemble(double density, double fovx, double fovy, int nx, int ny, double delta_t, double temp, 
                double len, double vis, int nrod,
                double link_len, vector<double *> pos_sets, double stretching, double bending, double frac_force, 
                string bc, double seed);
        
        void update_stretching();
        
        void update_bending();
        
};
  
class DLfilament_ensemble:
    public filament_ensemble<DLfilament>
{
    public:
        DLfilament_ensemble(double density, double fovx, double fovy, int nx, int ny, double delta_t, double temp, 
                double len, double vis, int nrod,
                double link_len, vector<double *> pos_sets, double stretching, double bending, double frac_force, 
                double bending_frac_force, string bc, double seed);
        
        void update_bending();
        
        void update_stretching();
        
        void set_bending_linear();
};

#endif
