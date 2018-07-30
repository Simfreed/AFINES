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
#include "filament.h"
//=====================================
//filament network class
class filament_ensemble
{
    public:
        
        filament_ensemble();

        filament_ensemble(int npolymer, int nbeads_min, int nbeads_max, double nbeads_prob,
                array<double,2> myfov, array<int,2> mynq, double delta_t, double temp,
                double rad, double vis, double spring_len, vector<array<double, 3> > pos_sets, double stretching, double ext, double bending, 
                double frac_force, string bc, double seed);

        filament_ensemble(double density, array<double,2> myfov, array<int, 2> mynq, double delta_t, double temp, 
                double len, double vis, int nbead,
                double spring_len, vector<array<double, 3> > pos_sets, double stretching, double ext, double bending, double frac_force, 
                string bc, double seed, double RMAX, double A);
        
        filament_ensemble(vector< vector<double> > beads, array<double,2> myfov, array<int,2> mynq, double delta_t, double temp,
                double vis, double spring_len, double stretching, double ext, double bending, double frac_force, string bc, double RMAX, double A); 
        
        ~filament_ensemble();
        
        void nlist_init();
        
        void nlist_init_serial();
        
        void quad_update();
        
        void quad_update_serial();
        
        void consolidate_quads();

        void update_quads_per_filament(int);

        void reset_n_springs(int);

        void update_dist_map(set<pair<double, array<int, 2>>>& t_map, const array<int, 2>& mquad, double x, double y);
        
        vector<filament *> * get_network();

        filament * get_filament(int index);

        set<pair<double, array<int,2>>> get_dist(double x, double y);
        
        set<pair<double, array<int,2>>> get_dist_all(double x, double y);
        
        array<double,2> get_direction(int fil, int spring);

        array<double,2> get_start(int fil, int spring);
        
        array<double,2> get_end(int fil, int spring);
        
        array<double,2> get_force(int fil, int bead);
        
        double get_int_direction(int fil, int spring, double xp, double yp);

        double get_xcm(int fil, int spring);
       
        double get_ycm(int fil, int spring);

        double get_llength(int fil, int spring);
       
        double get_bead_friction();
        
        double get_delrx();
        
        double get_stretching_energy();
        
        double get_bending_energy();

        double get_exv_energy();
        
        double get_kinetic_energy_vir(); 
        
        int get_nbeads();
        
        int get_nsprings();
        
        int get_nfilaments();

        vector<vector<double> > spring_spring_intersections(double cllen, double prob);

        void update_shear();
        
        void update_d_strain(double);
        
        void update_delrx(double);
        
        void update_stretching();
        
        void update_filament_stretching(int);
        
        void update_bending();
        
        void update_int_forces();

        void update_positions();

        void update_positions_range(int lo, int hi);
        
        void update_forces(int fil, int bead, double f2, double f3);

        void update_spring_forces(int f); 

        void update_spring_forces_from_quads(); 

        void update_force_between_filaments(double n1, double l1, double n2, double l2); 

        void update_excluded_volume(int f); 

        void write_beads(ofstream& fout);
        
        void write_springs(ofstream& fout);
        
        void write_thermo(ofstream& fout);
        
        void set_straight_filaments(bool is_straight);

        void set_y_thresh(double);
        
        void set_fene_dist_pct(double);
        
        void set_shear_rate(double);
        
        void set_shear_stop(double);

        void set_shear_dt(double);
        
        bool is_polymer_start(int f, int a);

        void set_fov(double x, double y);

        void set_nq(double x, double y);

        void set_visc(double v);

        vector<int> get_broken();

        void clear_broken();
        
        void print_filament_thermo();

        void print_network_thermo();

        void print_filament_lengths();
        
        void update();
        
        void update_energies();
        
        void turn_quads_off();
        
        void set_growing(double, double, double, double, int);

    protected:

        double t, dt, temperature, spring_rest_len, visc, min_time;
        double gamma, shear_stop, shear_dt, shear_speed, delrx;
        int nsprings_per_fil_max;
        
        bool straight_filaments = false, quad_off_flag;
        double pe_stretch, pe_bend, pe_exv, ke_vir;
        string BC;
        double rmax; 
        double kexv;  

        /* Here above is our excluded volume constant */

        array<double,2> fov, view;
        array<int, 2> nq, half_nq;
        vector<int> broken_filaments, empty_vector;
        
        vector< vector < vector< array<int, 2 > >* > * > springs_per_quad;
        vector< vector < int >* > n_springs_per_quad;
        
        vector<array<int, 2>* > all_quads;
        vector<filament *> network;
        unordered_set<array<int, 2>, boost::hash<array<int,2>>> fls;
};

#endif
