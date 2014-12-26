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
class filament_ensemble
{
    public:
        filament_ensemble();

        filament_ensemble(double density, double fovx, double fovy, int nx, int ny, double delta_t, double temp, 
                double len, double vis, int nrod,
                double link_len, std::vector<double *> pos_sets, double stretching, double bending, double seed);
        
        ~filament_ensemble();
        
        void quad_update_monomer(int i);
        
        void quad_update();

        std::vector<filament *> * get_network();

        std::map<std::vector<int>, double> get_dist(double x, double y);

        double* get_direction(int fil, int rod);

        double* get_intpoints(int fil, int rod, double xp, double yp);

        double get_int_direction(int fil, int rod, double xp, double yp);

        double get_xcm(int fil, int rod);
        
        double get_ycm(int fil, int rod);

        double get_angle(int fil, int rod);

        double get_alength(int fil, int rod);

        double* get_start(int fil, int rod);
        
        double* get_end(int fil, int rod);
        
        double* get_forces(int fil, int rod);

        void update(double t);

        void update_forces(int fil, int rod, double f1, double f2, double f3);

        void write(std::ofstream& fout);
        
        std::vector<double> get_angle_correlation(int polymer_index);

        std::map<int, std::vector<double> > get_all_angle_correlations();

        double get_fourier_mode(int n, int polymer_index);

        void connect_polymers(link_ensemble * links, double link_length, double link_stiffness, double bending_stiffness, std::string link_color);

        void update_polymer_bending(int polymer_index);

        void update_bending();
        
        void update_polymer_excluded_volume(int polymer_index);

        void update_excluded_volume();
        
        void add_monomer(filament * a, int n);

        void clear_filament_link_map();
        
        void set_straight_filaments(bool is_straight);

        void set_shear_rate(double);

        void update_polymer_shear(int polymer_index);

        void update_shear();
        
        bool is_polymer_start(int f, int r);

        void set_fov(double x, double y);

        void set_nq(double x, double y);

        void set_ld(double l);

        void set_visc(double v);

    private:
        double dt, temperature, fov[2], view[2], rho, ld, visc, gamma;
        
        double link_ld;
        
        int npolymer, nmon, nq[2];
       
        bool straight_filaments = false;

        std::vector<filament *> network;
        
        std::vector<int> empty_vector;
        
        std::map<int, std::map<int,std::vector< std::vector<int> > > > quad_fils;
        
        std::map<std::vector<int>,double> t_map;
};
#endif
