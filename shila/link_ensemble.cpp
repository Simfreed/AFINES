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
        fout<<links[i]->get_heads()[0]<<"\t"<<links[i]->get_heads()[1]<<"\t"
            <<links[i]->get_heads()[2]-links[i]->get_heads()[0]<<"\t"<<
            links[i]->get_heads()[3]-links[i]->get_heads()[1]<<"\t"<<
            links[i]->get_color()<<"\n";
    } 
}

void link_ensemble::add_link(Link * l)
{
    links.push_back(l);
}


void link_ensemble::delete_all()
{
    for (unsigned int i=0; i<links.size(); i++){
        delete links[i];
    }
}
