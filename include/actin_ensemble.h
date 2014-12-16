/*
 *  actin.cpp
 *  
 *
 *  Created by Shiladitya Banerjee on 9/3/13.
 *  Copyright 2013 University of Chicago. All rights reserved.
 *
 */

//=====================================
//include guard
#ifndef __ACTIN_ENSEMBLE_H_INCLUDED__
#define __ACTIN_ENSEMBLE_H_INCLUDED__

//=====================================
// forward declared dependencies
class Link;
class link_ensemble;

//=====================================
//included dependences
#include "string"
#include "vector"
#include "map"
#include "actin.h"

//=====================================
//actin network class
class actin_ensemble
{
    public:
        actin_ensemble();

        actin_ensemble(double density, double fovx, double fovy, int nx, int ny, double delta_t, double temp, 
                double len, double vis, int nmonomer,
                double link_len, std::vector<double *> pos_sets, double seed);
        
        ~actin_ensemble();
        
        void add_polymer(double startx, double starty, double theta, int index, int nmonomers);

        void quad_update_monomer(int i);
        
        void quad_update();

        std::vector<actin *> * get_network();

        std::map<int,double> get_dist(double x, double y);

        double* get_direction(int index);

        double* get_intpoints(int index, double xp, double yp);

        double get_int_direction(int index, double xp, double yp);

        double get_xcm(int index);
        
        double get_ycm(int index);

        double get_angle(int index);

        double get_alength(int index);

        double* get_start(int index);
        
        double* get_end(int index);
        
        double* get_forces(int index);

        void update(double t);

        void update_forces(int index, double f1, double f2, double f3);

        void write(std::ofstream& fout);
        
        std::vector<double> get_angle_correlation(int polymer_index);

        std::map<int, std::vector<double> > get_all_angle_correlations();

        double get_fourier_mode(int n, int polymer_index);

        void connect_polymers(link_ensemble * links, double link_length, double link_stiffness, double bending_stiffness, std::string link_color);

        void update_polymer_bending(int polymer_index);

        void update_bending();
        
        void update_polymer_excluded_volume(int polymer_index);

        void update_excluded_volume();
        
        void add_monomer(actin * a, int n);

        void clear_actin_link_map();
        
        void set_straight_filaments(bool is_straight);

        void set_shear_rate(double);

        void update_polymer_shear(int polymer_index);

        void update_shear();
        
        bool is_polymer_start(int i);

        void set_fov(double x, double y);

        void set_nq(double x, double y);

        void set_ld(double l);

        void set_visc(double v);

    private:
        double dt, temperature, fov[2], view[2], rho, ld, visc, gamma;
        
        double link_ld;
        
        int npolymer, nmon, nq[2];
       
        bool straight_filaments = false;

        std::vector<actin *> network;
        
        std::vector<int> empty_vector;
        
        std::map<int, std::map<int, Link * > > actin_link_map; //Maps {aindex0, aindex1} --> *Link
        
        std::map<int, std::vector<int> > mono_map; //Maps {filament index} --> {list of monomer indices}        
        
        std::map<int, std::map<int,std::vector<int> > > quad_fils;
        
        std::map<int,double> t_map;
};
#endif
