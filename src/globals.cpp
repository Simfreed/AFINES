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
#include <boost/range/irange.hpp>
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

array<double, 2> rij_periodic(double dx, double dy, double xbox, double ybox)
{
    //Using the minimum image convention
    //Allen and Tildesley, page 30
    double rxij = dx - xbox * round(dx / xbox);
    double ryij = dy - ybox * round(dy / ybox);
    return {rxij, ryij};
}

array<double, 2> rij_xperiodic(double dx, double dy, double xbox, double ybox)
{
    //Using the minimum image convention
    //Allen and Tildesley, page 30
    double rxij = dx - xbox * round(dx / xbox);
    double ryij = dy;
    return {rxij, ryij};
}

array<double, 2> rij_lees_edwards(double dx, double dy, double xbox, double ybox, double delrx)
{
    //Using the minimum image convention
    //Allen and Tildesley, page 247 
    double cory, rxij, ryij;
    cory = round(dy / ybox);
    rxij = dx   - cory * delrx;
    rxij = rxij - round(rxij / xbox) * xbox;
    ryij = dy   - cory * ybox;
    return {rxij, ryij};
}

double dist_bc(string bc, double dx, double dy, double xbox, double ybox, double delrx){
    
    array<double, 2> rij = rij_bc(bc, dx, dy, xbox, ybox, delrx);
    return hypot(rij[0], rij[1]);
}

array<double, 2> rij_bc(string bc, double dx, double dy, double xbox, double ybox, double delrx){
    
    if (bc == "PERIODIC")
        return rij_periodic(dx, dy, xbox, ybox);
    else if (bc == "XPERIODIC")
        return rij_xperiodic(dx, dy, xbox, ybox);
    else if (bc =="LEES-EDWARDS")
        return rij_lees_edwards(dx, dy, xbox, ybox, delrx);
    else
        return {dx, dy};

}

double dot_bc(string bc, double dx1, double dy1, double dx2, double dy2, double xbox, double ybox, double delrx)
{
    array<double, 2> 
        rij1 = rij_bc(bc, dx1, dy1, xbox, ybox, delrx), 
        rij2 = rij_bc(bc, dx2, dy2, xbox, ybox, delrx);
    
    return dot(rij1[0], rij1[1], rij2[0], rij2[1]);
}

double velocity(double vel0, double force, double fstall)
{
    double v;
    if (force>=fstall) {
        v = 0;
    }
    else if (force>0 && force<fstall) {
        v = vel0*(1-(fabs(force)/fstall));
    }
    else if (force<=0 && force>=-fstall) {
        v = vel0*(1+(fabs(force)/fstall));
    }
    else{
        v =  2*vel0;
    }
    return v;
}

array<double, 2> cm_bc(string bc, vector<double> xi, vector<double> yi, double xbox, double ybox, double delrx)
{
    if (bc == "PERIODIC" || bc == "LEES-EDWARDS")
        return {mean_periodic(xi, xbox) , mean_periodic(yi, ybox)};
    else
        return {mean(xi), mean(yi)};
}

double mean(vector<double> nums)
{
    double tot = 0;
    for (double n : nums) tot += n;
    return tot/((double) nums.size());
}

// Source https://en.wikipedia.org/wiki/Center_of_mass#Systems_with_periodic_boundary_conditions
double mean_periodic(vector<double> nums, double bnd)
{
    double theta, xitot=0, zetatot=0;
    for (double n : nums){
        theta = n*2*pi/bnd;
        xitot += cos(theta);
        zetatot += sin(theta);
    }
    double thetabar = atan2(zetatot/((double) nums.size()), xitot/((double) nums.size())) + pi; 
    return bnd*thetabar/(2*pi);
}

double cross(double ax, double ay, double bx, double by)
{
    return ax*by-bx*ay;
}

double dot(double x1, double y1, double x2, double y2)
{
    return x1*x2+y1*y2;
}

double dot(array<double, 2> v1, array<double, 2> v2)
{
    return v1[0]*v2[0]+v1[1]*v2[1];
}

double var(vector<double> vals)
{
    double m = mean(vals), sum = 0;
    for (unsigned int i = 0; i < vals.size(); i++){
        sum += (vals[i] - m)*(vals[i] - m);
    }
    return sum / vals.size();
}

double mode_var(vector<double> vals, double m)
{
    double sum = 0;
    for (unsigned int i = 0; i < vals.size(); i++){
        sum += (vals[i] - m)*(vals[i] - m);
    }
    return sum / vals.size();
}

vector<double> sum_vecs(vector<double> v1, vector<double> v2)
{
    vector<double> s;
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

vector<double *> vec2ptrvec(vector<double> v, int dim)
{
    vector<double *> out;
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

vector<double *> str2ptrvec(string pos_str, string pos_dlm, string coord_dlm)
{
    vector<string> posn, posns, coords;           
    vector<double *> out;
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

vector<array<double,3> > str2arrvec(string pos_str, string pos_dlm, string coord_dlm)
{
    vector<string> posn, posns, coords;           
    vector<array<double,3> > out;
    array<double,3> pos;

    boost::split(posns, pos_str, boost::is_any_of(pos_dlm));

    for(unsigned int i=0; i < posns.size(); i++){
        
        boost::split(coords, posns[i], boost::is_any_of(coord_dlm));
        
        for(unsigned int j=0; j < 3; j++) pos[j] = (double) atof(coords[j].data());
        
        out.push_back(pos);
    }

    return out;
}

vector<vector<double> > file2vecvec(string path, string delim)
{
    vector<vector<double> > out;
    string pos_str = "";
    vector<string> coords;
    vector<double> pos;
    
    ifstream pos_file;
    pos_file.open(path);
    
    while(getline(pos_file, pos_str))
    {
        boost::trim_right(pos_str);
        boost::split(coords, pos_str, boost::is_any_of(delim));
        
        for(unsigned int j=0; j < coords.size(); j++) 
            pos.push_back( (double) atof(coords[j].data()) );
        
        out.push_back(pos);

        pos.clear();
    }

    pos_file.close();
    
    return out;
}

void intarray_printer(array<int, 2> a)
{
    cout<<"\n{ " <<a[0]<<" , "<<a[1]<<" }";
}

template <typename T> int sgn(T val){
    return (T(0) < val) - (val < T(0));
}

vector<int> int_range(int lo, int hi)
{
    vector<int> out;
    for (int i = lo; i< hi; i++) out.push_back(i);
    return out;
}

vector<int> range_bc(string bc, double delrx, int topq, int lo, int hi)
{
    vector<int> out;
    if (lo <= hi)
        out = int_range(lo, hi);
    else if (bc == "PERIODIC" || bc == "LEES-EDWARDS"){
        vector<int> A = int_range(lo, topq), B = int_range(-topq, hi);
        out.reserve(A.size() + B.size());
        out.insert(out.end(), A.begin(), A.end());
        out.insert(out.end(), B.begin(), B.end());
    }
    else 
        out = vector<int>();

    return out;
}

array<double, 2> pos_bc(string bc, double delrx, double dt, array<double, 2> fov, array<double, 2> vel, array<double, 2> pos)
{
    double xnew = pos[0], ynew = pos[1];
        
    double xleft  = -fov[0] * 0.5;
    double xright =  fov[0] * 0.5;
    double yleft  = -fov[1] * 0.5;
    double yright =  fov[1] * 0.5;

    if(bc == "REFLECTIVE")
    {
        double local_shear = delrx * 2 * ynew / fov[1];
        xleft  += local_shear; //sheared simulation bounds
        xright += local_shear;
        if (xnew <= xleft || xnew >= xright) xnew -= 2*dt*vel[0];
        if (ynew <= yleft || ynew >= yright) ynew -= 2*dt*vel[1];

    }
    else if(bc == "INFINITE")
    {
        double local_shear = delrx * 2 * ynew / fov[1];
        xleft  += local_shear; //sheared simulation bounds
        xright += local_shear;
        if      (xnew <= xleft)  xnew = xleft;
        else if (xnew >= xright) xnew = xright;
        if      (ynew <= yleft)  ynew = yleft;
        else if (ynew >= yright) ynew = yright;

    }
    else if(bc == "XPERIODIC")
    {
        if      (xnew < xleft)  xnew += fov[0];
        else if (xnew > xright) xnew -= fov[0];

        if      (ynew < yleft)  ynew = yleft;
        else if (ynew > yright) ynew = yright;

    }
    else if(bc == "PERIODIC")
    {
        xnew = xnew - fov[0] * round(xnew / fov[0]);
        ynew = ynew - fov[1] * round(ynew / fov[1]);
    }
    else if(bc == "LEES-EDWARDS")
    {
        double cory = round(ynew/fov[1]);
        xnew = xnew - delrx  * cory;
        xnew = xnew - fov[0] * round(xnew / fov[0]);
        ynew = ynew - fov[1] * cory;
    }
    
    return {xnew, ynew};

}


// Method to sort a map by value; source, for more general formulation:
// http://stackoverflow.com/questions/5056645/sorting-stdmap-using-value/5056797#5056797
pair<double, array<int, 2> > flip_pair(const pair<array<int, 2>, double> &p)
{
        return std::pair<double,array<int,2> >(p.second, p.first);
}

multimap<double, array<int, 2> > flip_map(const map<array<int, 2>, double> &src)
{
    multimap<double,array<int,2> > dst;
    std::transform(src.begin(), src.end(), std::inserter(dst, dst.begin()), 
            flip_pair);
    return dst;
}


template int sgn<int>(int);
template int sgn<double>(double);
template int sgn<float>(float);
