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
int event(double rate, double timestep);
double dis_points(double x1, double y1, double x2, double y2);
double velocity(double vel0, double force, double fstall);
double cross(double ax, double ay, double bx, double by);
double dot(double x1, double y1, double x2, double y2);
double mean(vector<double> vals);
double var(vector<double> vals);
double mode_var(vector<double> vals, double m);
bool close(double e, double a, double r);
vector<double> sum_vecs(vector<double> v1, vector<double> v2);
vector<double *> vec2ptrvec(vector<double>, int dim);
vector<double *> str2ptrvec(string, string, string);
vector<array<double,3> > str2arrvec(string, string, string);
vector<vector<double> > file2vecvec(string path, string delim);

map<array<int, 2>, double> transpose(map<array<int, 2>, double> mat);
map<array<int, 2>, double> invert_block_diagonal(map<array<int, 2>, double> mat);
void intarray_printer(array<int,2> a);
#endif
