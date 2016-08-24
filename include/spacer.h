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
                array<double, 2> myfov, double delta_t, double v0, double temp, double stiffness, double max_ext_ratio, 
                double ron, double roff, double rend, 
                double fstall, double fbreak, double bindEng,
                double vis, string BC);
        
        spacer(array<double, 4> pos, double mlen, filament_ensemble* network, 
                array<int, 2> mystate, array<int, 2> myfindex, array<int, 2> myrindex,
                array<double, 2> myfov, double delta_t, double v0, double temp, double stiffness, double max_ext_ratio,
                double ron, double roff, double rend, 
                double fstall, double fbreak, double bindEng,
                double vis, string BC);

        ~spacer();

        void set_bending(double, double);

        double get_kb();
        
        double get_th0();

        void update_force();
        
        void identify();
        
        void update_bending(int);

    private:

        double kb, th0; 
        
};

#endif
