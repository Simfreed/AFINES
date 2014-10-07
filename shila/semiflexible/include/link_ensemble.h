/*
 *  link_ensemble.h
 *
 *  Created by Simon Freedman on 9/10/2014
 *  Copyright 2014 University of Chicago. All rights reserved.
 *
 */

//=====================================
//include guard
#ifndef __LINK_ENSEMBLE_H_INCLUDED__
#define __LINK_ENSEMBLE_H_INCLUDED__

//=====================================
// forward declared dependencies

//=====================================
//included dependences
#include "vector"
#include "Link.h"

//=====================================
//link ensemble class
class link_ensemble
{
    public:

        link_ensemble();
        ~link_ensemble();

        void link_walk();

        void link_write(std::ofstream& fout);

        void add_link(Link * l);
        
        void clear();
        
        int size();

        Link * at(int index);
    protected:
        std::vector<Link *> links;
};
#endif
