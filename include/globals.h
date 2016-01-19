/*
 *  globals.h
 *  
 *  Created by Simon Freedman on 9/12/2014
 *
 *  Created by Shiladitya Banerjee on 9/3/13.
 *  Copyright 2013 University of Chicago. All rights reserved.
 *
 */

//=====================================
//include guard
#ifndef __GLOBALS_H_INCLUDED__
#define __GLOBALS_H_INCLUDED__

//=====================================
// forward declared dependencies

//=====================================
//included dependences
#include <iostream> //std::cout
#include <boost/algorithm/string.hpp>
#include <boost/range/irange.hpp>
#include <boost/range/join.hpp>
#include <boost/optional.hpp>
//#include "iomanip"
#include <math.h>
#include "fstream"
#include "string"
#include "stdint.h"
#include <stdlib.h>
#include "vector"
#include "array"
#include <stdio.h>
#include <array>
#include <map>
#include <algorithm> //std::for_each
#include <unordered_set>
#include <limits>
#include <cstddef>

using namespace std;

/* distances in microns, time in seconds, forces in pN * 
 * --> Temp in pN-um                                   */
//const double pi, eps, dt, temperature;
const double pi = 3.14159265358979323;
const double maxSmallAngle = pi/12; //Small angles DEFINED as such that sin(t) = t to 2 SigFigs
const double eps = 1e-10;
const double temperature = 0.004;
const double actin_mass_density = 2.6e-14; //miligram / micron
/*generic functions to be used below*/

double rng(double start, double end);
int pr(int num);
double rng_exp(double mean);
double rng_n(double mean, double var);
bool event(double prob);
int event(double rate, double timestep);

array<double, 2> rij_periodic(double dx, double dy, double xbox, double ybox);
array<double, 2> rij_lees_edwards(double dx, double dy, double xbox, double ybox, double shear_dist);
array<double, 2> rij_bc(string bc, double dx, double dy, double xbox, double ybox, double shear_dist);

vector<int> range_bc(string bc, double delrx, int topq, int low, int high);
vector<int> int_range(int lo, int hi);

double mean_periodic(vector<double> nums, double bnd);
double mean(vector<double> nums);
array<double, 2> cm_bc(string bc, vector<double> xi, vector<double> yi, double xbox, double ybox, double shear_dist);

double dist_bc(string bc, double dx, double dy, double xbox, double ybox, double shear_dist);
double dot_bc(string bc, double dx1, double dy1, double dx2, double dy2, double xbox, double ybox, double shear_dist);
array<double, 2> pos_bc(string bc, double delrx, double dt, array<double, 2> fov, array<double, 2> vel, array<double, 2> pos);

double velocity(double vel0, double force, double fstall);
double cross(double ax, double ay, double bx, double by);
double dot(double x1, double y1, double x2, double y2);
double dot(array<double, 2> v1, array<double, 2> v2);

double var(vector<double> vals);
double mode_var(vector<double> vals, double m);
bool close(double e, double a, double r);
bool are_same(double a, double b);
vector<double> sum_vecs(vector<double> v1, vector<double> v2);
vector<double *> vec2ptrvec(vector<double>, int dim);
vector<double *> str2ptrvec(string, string, string);
vector<array<double,3> > str2arrvec(string, string, string);
vector<vector<double> > file2vecvec(string path, string delim);

template <typename T> int sgn(T val);

pair<double, array<int, 2> > flip_pair(const pair<array<int, 2>, double> &p);
multimap<double, array<int, 2> > flip_map(const map<array<int, 2>, double> &p);
string print_pair(string name, array<double, 2> p);

//template <typename A, typename B> pair<B,A> flip_pair(const pair<A,B> &p);
//template <typename A, typename B> multimap<B,A> flip_map(const map<A,B> &src);

map<array<int, 2>, double> transpose(map<array<int, 2>, double> mat);
map<array<int, 2>, double> invert_block_diagonal(map<array<int, 2>, double> mat);
void intarray_printer(array<int,2> a);

boost::optional<array<double, 2> > seg_seg_intersection(array<double, 2>, array<double, 2>, array<double, 2>, array<double, 2>);
boost::optional<array<double, 2> > seg_seg_intersection_bc(string, double, array<double, 2>, array<double, 2>, array<double, 2>, array<double, 2>, array<double, 2>);
#endif
