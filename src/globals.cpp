/*
 *  globals.cpp
 *  
 *
 *  Adapted by Simon Freedman on 9/12/14
 *  Created by Shiladitya Banerjee on 9/3/13.
 *  Copyright 2013 University of Chicago. All rights reserved.
 *
 */

#include "globals.h"

/* distances in microns, time in seconds, forces in pN */

/*generic functions to be used below*/
double rng(double start, double end)
{
	return start+(end-start)*((double)rand()/(RAND_MAX));
}

int pr(int num)
{
	if (num==0) {
		return 1;	
	}
	else {
		return 0;
	}
}

double rng_exp(double mean)
{
    double u;
    u=rand() / (RAND_MAX + 1.);
    return  -mean*log(u);
}


double rng_n(double mean, double var)
{
    double U = 0, V, Z;
    
    while(U==0) U = rng(0,1);
    V = rng(0,1);
    Z = sqrt(-2 * log(U)) * sin(2 * pi * V);
    
    return mean+var*Z;
}

int event(double rate, double timestep)
{
    if (rng(0,1.0)<rate*timestep) {
        return 1;
    }
    else
        return 0;
}

double dis_points(double x1, double y1, double x2, double y2)
{
    double dis=sqrt(pow(x2-x1,2)+pow(y2-y1,2));
    return dis;
}

double velocity(double vel0, double force, double fstall)
{
    if (force>=fstall) {
        return 0;
    }
    else if (force>0 && force<fstall) {
        return vel0*(1-(fabs(force)/fstall));
    }
    else if (force<=0 && force>=-fstall) {
        return vel0*(1+(fabs(force)/fstall));
    }
    else{
        return 2*vel0;
    }
}

double cross(double ax, double ay, double bx, double by)
{
    return ax*by-bx*ay;
}

double dot(double x1, double y1, double x2, double y2)
{
    return x1*x2+y1*y2;
}

double mean(std::vector<double> vals)
{
    double sum = 0;
    for (unsigned int i = 0; i < vals.size(); i++){
        sum+= vals[i];
    }
    return sum / vals.size();

}

double var(std::vector<double> vals)
{
    double m = mean(vals), sum = 0;
    for (unsigned int i = 0; i < vals.size(); i++){
        sum += (vals[i] - m)*(vals[i] - m);
    }
    return sum / vals.size();
}

double mode_var(std::vector<double> vals, double m)
{
    double sum = 0;
    for (unsigned int i = 0; i < vals.size(); i++){
        sum += (vals[i] - m)*(vals[i] - m);
    }
    return sum / vals.size();
}

std::vector<double> sum_vecs(std::vector<double> v1, std::vector<double> v2)
{
    std::vector<double> s;
    if (v1.empty())
        s = v2;
    else if( v2.empty())
        s = v1;
    else if (v1.size() != v2.size())
        return s;
    else{
        for (unsigned int i = 0; i < v1.size(); i++){
            s.push_back(v1[i] + v2[i]);
        }
    }
    return s;
}

bool close(double actual, double expected, double err)
{
    if (expected == 0){
        return fabs(expected-actual) < err;
    }
    else{
        return fabs(expected-actual)/expected < err;
    }
}

/* Takes a vector formatted 
 * [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
 * And converts it to a vector formatted
 * [{1, 2, 3, 4}, {5, 6, 7, 8}, {9, 10, 11, 12}] if dim = 4 
 * Or
 * [{1, 2, 3}, {4, 5, 6}, {7, 8, 9}, {10, 11, 12}] if dim = 3
 */

std::vector<double *> vec2ptrvec(std::vector<double> v, int dim)
{
    std::vector<double *> out;
    double * pos;
    for (unsigned int i = 0; i < v.size(); i+=dim)
    {
        pos = new double[dim];
        for (int j = 0; j < dim; j++)
        {
           pos[j] = v[i+j];
        }
        out.push_back(pos);
    }
    return out;
}

/* Takes a string formatted
 * 1,2,3;4,5,6;7,8,9;10,11,12
 * and converts it into a vector of pointers:
 * [{1, 2, 3}, {4, 5, 6}, {7, 8, 9}, {10, 11, 12}] 
 */

std::vector<double *> str2ptrvec(std::string pos_str, std::string pos_dlm, std::string coord_dlm)
{
    std::vector<std::string> posn, posns, coords;           
    std::vector<double *> out;
    double * pos;

    boost::split(posns, pos_str, boost::is_any_of(pos_dlm));

    for(unsigned int i=0; i < posns.size(); i++){
        
        boost::split(coords, posns[i], boost::is_any_of(coord_dlm));
        pos = new double[coords.size()];
        
        for(unsigned int j=0; j < coords.size(); j++){
            pos[j] = (double) atof(coords[j].data());
        }
        out.push_back(pos);
    }

    return out;
}
