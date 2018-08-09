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
#include <boost/functional/hash.hpp>
#define BOOST_NO_CXX11_SCOPED_ENUMS
#include <boost/filesystem.hpp>
#undef BOOST_NO_CXX11_SCOPED_ENUMS

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
#include <unordered_map>
#include <algorithm> //std::for_each
#include <unordered_set>
#include <set>
#include <limits>
#include <cstddef>
// #include <omp.h>
#include <random>

using namespace std;
namespace fs = boost::filesystem;

/* distances in microns, time in seconds, forces in pN * 
 * --> Temp in pN-um                                   */
//const double pi, eps, dt, temperature;
const double pi = 3.14159265358979323;
const double maxSmallAngle = 0.001; //pi/12.0; //Small angles DEFINED as such that sin(t) = t to 2 SigFigs
const double eps = 1e-4;
const double infty = 1e10;
const double actin_mass_density = 2.6e-14; //miligram / micron
/*generic functions to be used below*/
void set_seed(int s);
double rng_u();
double rng(double start, double end);
int pr(int num);
double rng_exp(double mean);
double rng_n(); //default parameters --> mu = 0, sig = 1
double rng_n(double mean, double var);
bool event(double prob);
int event(double rate, double timestep);

array<double, 2> rij_periodic(double dx, double dy, double xbox, double ybox);
array<double, 2> rij_lees_edwards(double dx, double dy, double xbox, double ybox, double shear_dist);
array<double, 2> rij_bc(string bc, double dx, double dy, double xbox, double ybox, double shear_dist);

vector<int> range_bc(string bc, double delrx, int botq, int topq, int low, int high);
vector<int> range_bc(string bc, double delrx, int botq, int topq, int low, int high, int di);
vector<int> int_range(int lo, int hi);
vector<int> int_range(int lo, int hi, int di);

double mean_periodic(const vector<double>& nums, double bnd);
double mean(const vector<double>& nums);
array<double, 2> cm_bc(string bc, const vector<double>& xi, const vector<double>& yi, double xbox, double ybox, double shear_dist);

double dist_bc(string bc, double dx, double dy, double xbox, double ybox, double shear_dist);
double dot_bc(string bc, double dx1, double dy1, double dx2, double dy2, double xbox, double ybox, double shear_dist);
array<double, 2> pos_bc(string bc, double delrx, double dt, const array<double, 2>& fov, const array<double, 2>& vel, const array<double, 2>& pos);
array<int, 2> coord2quad(const array<double, 2>& fov, const array<int, 2>& nq, const array<double, 2>& pos);
int coord2quad_floor(double fov, int nq, double pos);
int coord2quad_ceil(double fov, int nq, double pos);
int coord2quad(double fov, int nq, double pos);

double my_velocity(double vel0, double force, double fstall);
double cross(double ax, double ay, double bx, double by);
double dot(double x1, double y1, double x2, double y2);
double cross(const array<double, 2>& v1, const array<double, 2>& v2);
double dot(const array<double, 2>& v1, const array<double, 2>& v2);

double angBC(double ang, double max);

double var(const vector<double>& vals);
double mode_var(const vector<double>& vals, double m);
bool close(double e, double a, double r);
bool are_same(double a, double b);
vector<double> sum_vecs(const vector<double>& v1, const vector<double>& v2);
vector<double *> vec2ptrvec(const vector<double>&, int dim);
vector<double *> str2ptrvec(string, string, string);
vector<array<double,3> > str2arrvec(string, string, string);
vector<vector<double> > file2vecvec(string path, string delim);
vector<vector<double> > traj2vecvec(string path, string delim, double tf);
double last_full_timestep(string dirpath);
void write_first_nlines(string path, int n);
void write_first_ntsteps(string path, int n);
void write_first_tsteps(string path, double tstop);

template <typename T> int sgn(T val);
int mysgn(double);

pair<double, array<int, 2> > flip_pair(const pair<array<int, 2>, double> &p);
multimap<double, array<int, 2> > flip_map(const std::unordered_map<array<int, 2>, double, boost::hash<array<int,2>>> &p);
string print_pair(string name, const array<double, 2>& p);

//template <typename A, typename B> pair<B,A> flip_pair(const pair<A,B> &p);
//template <typename A, typename B> multimap<B,A> flip_map(const map<A,B> &src);

map<array<int, 2>, double> transpose(map<array<int, 2>, double> mat);
map<array<int, 2>, double> invert_block_diagonal(map<array<int, 2>, double> mat);
void intarray_printer(array<int,2> a);

boost::optional<array<double, 2> > seg_seg_intersection(const array<double, 2>&, const array<double, 2>&, const array<double, 2>&, const array<double, 2>&);
boost::optional<array<double, 2> > seg_seg_intersection_bc(string, double, const array<double, 2>&, const array<double, 2>&, const array<double, 2>&, const array<double, 2>&, const array<double, 2>&);
std::string quads_error_message(std::string, vector<array<int, 2> >, vector<array<int, 2> > );
#endif



