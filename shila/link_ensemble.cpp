/*
 *  Link.cpp
 *  
 *
 *  Created by Simon Freedman on 9/10/2014
 *  Copyright 2014 University of Chicago. All rights reserved.
 *
 */

#include "link_ensemble.h"
#include "globals.h"

//link ensemble definitions 

link_ensemble::link_ensemble(){};
link_ensemble::~link_ensemble(){};

void link_ensemble::link_walk()
{
    for (unsigned int i=0; i<links.size(); i++) {
        links[i]->step();
        links[i]->actin_update();
    }
}

void link_ensemble::link_write(std::ofstream& fout)
{
    for (unsigned int i=0; i<links.size(); i++) {
        fout<<links[i]->to_string();
    } 
}

void link_ensemble::add_link(Link * l)
{
    links.push_back(l);
}

int link_ensemble::size()
{
    return links.size();
}

void link_ensemble::delete_all()
{
    for (unsigned int i=0; i<links.size(); i++){
        delete links[i];
    }
}
