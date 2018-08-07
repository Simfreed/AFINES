/*
 * spacer.h
 *  
 *
 *  Created by Simon Freedmane on 8/8/16.
 *  Copyright 2016 University of Chicago. All rights reserved.
 *
 */

//=====================================
//include guard
#ifndef __SPACER_H_INCLUDED__
#define __SPACER_H_INCLUDED__

//=====================================
// forward declared dependencies

//=====================================
//included dependences
#include "motor.h"

//spacer class

class spacer : public motor {
    public:

        spacer(array<double, 3> pos, double mlen, filament_ensemble* network, 
                array<int, 2> mystate, array<int, 2> myfindex, array<int, 2> myrindex,
                array<double, 2> myfov, double delta_t, double temp, double v0, double stiffness, double max_ext_ratio, 
                double ron, double roff, double rend, 
                double fstall, double rcut,
                double vis, string BC);
        
        spacer(array<double, 4> pos, double mlen, filament_ensemble* network, 
                array<int, 2> mystate, array<int, 2> myfindex, array<int, 2> myrindex,
                array<double, 2> myfov, double delta_t, double temp, double v0, double stiffness, double max_ext_ratio,
                double ron, double roff, double rend, 
                double fstall, double rcut,
                double vis, string BC);

        virtual ~spacer();

        void set_bending(double, double);

        double get_kb();
        
        double get_th0();

        virtual void update_force();
        
        void identify();
        
        int get_further_end(int, int, int);

        array<double, 2> disp_from_bead(int, int, int);

        void update_bending(int);
        
        array<array<double, 2>,2> get_b_force();

        virtual void brownian_relax(int hd);

        void filament_update();

        virtual double metropolis_prob(int hd, array<int, 2> flidx, array<double, 2> newpos, double maxRate);
        
        virtual bool allowed_bind(int hd, array<int, 2> flidx);
        
        array<double,2> get_bending_energy();

    private:

        double kb, th0; 
        array<array<double, 2>,2 > b_force;        
        array<double, 2> b_eng;
};

#endif
