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
    u=rand() / (RAND_MAX + 1.);;
    return  -mean*log(u);
}


double rng_n(double mean, double var)
{
    double U, V, Z;
    int phase = 0;
    
    if(phase == 0) {
        U = (rand() + 1.) / (RAND_MAX + 2.);
        V = rand() / (RAND_MAX + 1.);
        Z = sqrt(-2 * log(U)) * sin(2 * pi * V);
    } else
        Z = sqrt(-2 * log(U)) * cos(2 * pi * V);

    phase = 1 - phase;

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
    for (int i = 0; i < vals.size(); i++){
        sum+= vals[i];
    }
    return sum / vals.size();

}

double var(std::vector<double> vals)
{
    double m = mean(vals), sum = 0;
    for (int i = 0; i < vals.size(); i++){
        sum += (vals[i] - m)*(vals[i] - m);
    }
    return sum / vals.size();
}

double mode_var(std::vector<double> vals, double m)
{
    double sum = 0;
    for (int i = 0; i < vals.size(); i++){
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
        for (int i = 0; i < v1.size(); i++){
            s.push_back(v1[i] + v2[i]);
        }
    }
    return s;
}
