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

//=====================================
//actin filament class
class filament
{
    public:

        filament(double x0, double y0, double phi0, int nrod, double fovx, double fovy, int nx, int ny, 
                double vis, double deltat, double temp, bool isStraight,
                double rodLength, double linkLength, double stretching, double bending, double fracture); 

        filament(std::vector<actin *> rodvec, double linkLength, double stretching_stiffness, double bending_stiffness, 
                double deltat, double temp, double fracture);
       
        filament();
        
        ~filament();
    
        void set_shear(double g);

        void update(double t);
        
        void update_bending();
        
        std::vector<filament *> update_stretching();

        void update_shear(); 
        
        actin * get_rod(int i);
        
        Link * get_link(int i);

        std::vector<std::vector<std::vector<int> > > get_quadrants();
        
        std::string write_rods();
        
        std::string write_links();
        
        std::string to_string();
        
        bool operator==(const filament& that);    
    
        std::vector<actin *> get_rods(unsigned int first, unsigned int last);
        
        std::vector<filament *> fracture(int node);
        
        void update_forces(int index, double f1, double f2, double f3);

    private:
        
        double fov[2], gamma, temperature, dt, fracture_force;
        
        int nq[2];

        std::vector<std::vector<int> > quads_filled; //vector of two vectors(x and y quadrants) of integers
        
        std::vector<actin *> rods;
        
        std::vector<Link *> lks;
       
        std::vector<int> tmp;
};

#endif

