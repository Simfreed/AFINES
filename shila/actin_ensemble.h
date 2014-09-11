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
        actin_ensemble(double density, double fovx, double fovy, int nx, int ny, double len, double vis, int nmonomer,
                double link_len);
        
        ~actin_ensemble();

        void quad_update();

        std::vector<actin>* get_network();

        std::map<int,double> get_dist(double x, double y);

        double* get_direction(int index);

        double* get_intpoints(int index, double xp, double yp);

        double get_int_direction(int index, double xp, double yp);

        double* get_position(int index);

        double get_alength(int index);

        double* get_ends(int index);

        void update();

        void update_forces(int index, double f1, double f2, double f3);

        void write(std::ofstream& fout);
        
        std::vector<double> get_angle_correlation(int polymer_index);

        std::map<int, std::vector<double> > get_all_angle_correlations();

        double get_fourier_mode(int n, int polymer_index);

        void connect_polymers(link_ensemble * links, double link_length, double link_stiffness, std::string link_color);
        
        void connect_polymers(link_ensemble * links, 
                double link_length, double link_stiffness, std::string link_color,
                double b_link_stiffness, std::string b_link_color);
    
    private:
        double fov[2], rho, ld, xnew, ynew, phinew, vpar, omega, vperp, vx, vy, alength, view, a_ends[4], av_vel, visc;
        
        double link_ld, lx, ly, ltheta;
        
        int npolymer, nmonomer_max, nmonomer_min, nq[2], xn, yn, qxcm, qycm;
        
        std::vector<actin> network;
        
        std::vector<int> empty_vector;
        
        std::map<int, std::vector<double> > link_map; //Maps {monomer} --> {link_xcm, link_ycm, link_angle} 
        
        std::map<int, std::vector<int> > mono_map; //Maps {filament index} --> {list of monomer indices}        
        
        std::map<int, std::map<int,std::vector<int> > > quad_fils;
        
        std::map<int,double> t_map;
};
#endif
