/*
 *  link.h
 *
 *  Created by Simon Freedman on 9/10/2014
 *  Copyright 2014 University of Chicago. All rights reserved.
 *
 */

//=====================================
//include guard
#ifndef __LINK_H_INCLUDED__
#define __LINK_H_INCLUDED__


//=====================================
// forward declared dependencies
class filament;

//=====================================
//included dependences
#include "string"
#include "vector"
#include "globals.h"
#include "motor.h"
//=====================================
//Link class
class Link
{
    public:
        Link();
        
        Link(double len, double stiffness, double max_ext, filament* f, array<int, 2> aindex, array<double, 2> fov, array<int, 2> nq);
        
        virtual ~Link();

        array<double, 2> get_hx();
        
        array<double, 2> get_hy();
        
        double get_kl();
        
        double get_length();
       
        double get_angle();

        double get_l0();
        
        double get_fene_ext();
        
        double get_stretching_energy();
        
        double get_xcm();
        
        double get_ycm();
        
        string to_string();
        
        string write(string bc, double shear_dist);
        
        void step();
        
        void step(string bc, double shear_dist);
        
        void filament_update();
        
        bool operator==(const Link& that);    
        
        bool is_similar(const Link& that);    

        void update_force(string bc, double shear_dist);
        
        void update_force_fraenkel_fene(string bc, double shear_dist);
        
        double get_stretching_energy_fene(string bc, double shear_dist);
        
        void update_force_marko_siggia(string bc, double shear_dist, double kToverA);

        array<double,2> get_force();

        void set_aindex(array<int,2> idx);
        
        void set_l0(double myl0);
        
        double get_distance(string bc, double shear_dist, double xp, double yp);

        double get_int_angle(double xp, double yp);
        
//      array<double,2> get_intpoint(string bc, double shear_dist, double xp, double yp);
        array<double,2> get_intpoint();

        void calc_intpoint(string bc, double shear_dist, double xp, double yp);
        
        vector<array<int,2> > get_quadrants();
       
        void quad_update(string bc, double shear_dist);
        
        array<double, 2> get_direction();
        
        array<double, 2> get_disp();
        
        array<double, 2> get_neg_disp();
        
        double get_max_ext();
        
        void inc_aindex();

        void set_mots(map<int, int>* mymots);

        void set_xlinks(map<int, int>* myxlinks);
        
        map<int, int> * get_xlinks();
    
        void update_length();
        
        // stuff for growing
        void add_mot(motor * mot, int hd);

        void remove_mot(motor * mot);

        int get_n_mots();

        motor * get_mot(int i);
        
        int get_mot_hd(int i);
        
        array<int,2> get_aindex();
        
        map<motor *, int> & get_mots();
    
        void remove_all_mots();
        
        void set_fov(double, double);

        void set_nq(int, int);

    protected:

        double xcm, ycm, phi, l0, kl, max_ext, eps_ext, llen;//, force;
       
        array<double,2> fov, hx, hy;
        array<double, 2> disp, force, intpoint;

        array<int, 2> nq, half_nq, aindex;
         
        filament *fil;
        
        vector< array<int,2> > quad; //vector of two vectors(x and y quadrants) of integers

        map<motor *, int> mots;
};
#endif
