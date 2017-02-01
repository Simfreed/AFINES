#!/usr/bin/env python

import numpy as np
from scipy.interpolate import Rbf
from scipy.stats import binned_statistic_2d
import os
import argparse

def calc_div(x,y, xk, yk, wkx, wky, eps):
    return -2.0/eps**2*np.sum((wkx*(x-xk)+wky*(y-yk))*np.exp((-1.0/eps**2)*((x-xk)**2+(y-yk)**2)))

parser = argparse.ArgumentParser()

parser.add_argument("directory", type=str, help="directory of txt_stack/actins.txt file")
parser.add_argument("--parts", type=str, help="particles to track the velocity of", default = 'actins')
parser.add_argument("--fov", type=int, nargs=2, help="simulation cell size", default=[50,50])
parser.add_argument("--rbfunc", type=str, help="radial basis function", default='gaussian')
parser.add_argument("--rbeps", type=float, help="epsilon for radial basis function (e.g., R0 if rbfunc='gaussian')",default=5)
parser.add_argument("--dr", type=float, nargs=2, help="mesh grid size", default=[1,1])
parser.add_argument("--dt", type=int, help="number of timesteps to use when calculating velocity", default=10)
parser.add_argument("--padding", type=float, help="padding in x and y directions in fraction of cell", default=1)
parser.add_argument("--nbins", type=int, help="number of bins to use in smoothing data", default=10);
parser.add_argument("--minpts", type=int, help="min # of points to use when calculating means", default=10);
parser.add_argument("--div_lsc", type=float, help="size of a bin to calculating total divergence", default=10);

def thresh_mean(myarr):
    if len(myarr) < args.minpts:
        return np.nan
    else:
        return np.mean(myarr)

args = parser.parse_args()

fovarr = np.array(args.fov);
xsz, ysz = args.fov

fname = '{0}/txt_stack/{1}.txt'.format(args.directory, args.parts)

dat    = open(fname, 'r')
header = dat.readline()
nparts = int(header.rstrip().split('=')[-1])
dat.close()

dat   = np.loadtxt(fname, comments='t', usecols=(0,1))
nts   = len(dat)/float(nparts)
xyt   = np.array(np.split(dat[0:(nparts*np.floor(nts)).astype(int)], np.floor(nts)))
drt   = xyt[args.dt:]-xyt[:-args.dt]
drpt  = (drt-np.around(drt/fovarr)*fovarr)

rxyt  = xyt + np.array([ xsz,    0])
lxyt  = xyt + np.array([-xsz,    0])
uxyt  = xyt + np.array([   0,  ysz])
dxyt  = xyt + np.array([   0, -ysz])
urxyt = xyt + np.array([ xsz,  ysz])
ulxyt = xyt + np.array([-xsz,  ysz])
drxyt = xyt + np.array([ xsz, -ysz])
dlxyt = xyt + np.array([-xsz, -ysz])

xyall = np.concatenate([xyt, rxyt, lxyt, uxyt, dxyt, urxyt, ulxyt, drxyt, dlxyt], axis=1)
uvall = np.concatenate([drpt/args.dt]*9, axis=1)

f      = args.padding+0.5
edgesx = np.linspace(-f*xsz, f*xsz, args.nbins+1)
edgesy = np.linspace(-f*ysz, f*ysz, args.nbins+1)
print 'imported data'
mxt = np.stack([binned_statistic_2d(xyall[t,:,0], xyall[t,:,1], xyall[t,:,0], statistic=thresh_mean, bins=[edgesx, edgesy])[0].flatten() 
    for t in range(len(drpt))])
myt = np.stack([binned_statistic_2d(xyall[t,:,0], xyall[t,:,1], xyall[t,:,1], statistic=thresh_mean, bins=[edgesx, edgesy])[0].flatten() 
    for t in range(len(drpt))])
mut = np.stack([binned_statistic_2d(xyall[t,:,0], xyall[t,:,1], uvall[t,:,0], statistic=thresh_mean, bins=[edgesx, edgesy])[0].flatten() 
    for t in range(len(drpt))])
mvt = np.stack([binned_statistic_2d(xyall[t,:,0], xyall[t,:,1], uvall[t,:,1], statistic=thresh_mean, bins=[edgesx, edgesy])[0].flatten() 
    for t in range(len(drpt))])

print 'calculated bin statistic'

mxyt = np.stack([mxt, myt], axis=2) 
muvt = np.stack([mut, mvt], axis=2)

xi = np.arange(-0.5*xsz, 0.5*xsz + args.dr[0], args.dr[0])
yi = np.arange(-0.5*ysz, 0.5*ysz + args.dr[1], args.dr[1])

xmid = 0.5*(xi[1:]+xi[:-1])
ymid = 0.5*(yi[1:]+yi[:-1])

partname = args.parts if args.parts != 'actins' else '';
divxydir='{0}/analysis/xyuvd_nbins{1}_lsc{2}_dt{3}_thresh{4}{5}'.format(args.directory,args.nbins,args.div_lsc,args.dt,args.minpts,partname) 
weightsdir='{0}/analysis/vfield_weights_nbins{1}_lsc{2}_dt{3}_thresh{4}{5}'.format(args.directory,args.nbins,args.div_lsc,args.dt,args.minpts,partname) 

div_bins = []
if not os.path.exists(divxydir):
    os.makedirs(divxydir)
if not os.path.exists(weightsdir):
    os.makedirs(weightsdir)

f      = 0.5
edgesx = np.arange(-f*xsz, f*xsz + args.div_lsc, args.div_lsc)
edgesy = np.arange(-f*ysz, f*ysz + args.div_lsc, args.div_lsc)

for t in range(len(mxyt)):

    idxs = np.logical_and( np.isfinite( mxyt[t,:,0] ) , np.isfinite(mxyt[t,:,1]) )
    
    if len(mxyt[t,idxs]) > 0:
        
        rbfu = Rbf(mxyt[t,idxs,0], mxyt[t,idxs,1], muvt[t,idxs,0], function=args.rbfunc, epsilon=args.rbeps)
        rbfv = Rbf(mxyt[t,idxs,0], mxyt[t,idxs,1], muvt[t,idxs,1], function=args.rbfunc, epsilon=args.rbeps)

        wu = rbfu.nodes
        wv = rbfv.nodes

        xywuwve=np.stack([mxyt[t,idxs,0].flatten(), mxyt[t,idxs,1].flatten(), wu.flatten(), wv.flatten(),
            args.rbeps*np.ones(wu.shape[0])], axis=1)
        np.savetxt('{0}/t{1}.txt'.format(weightsdir, t), xywuwve, '%10.5f')

        print 't = {0}'.format(t)
    else:
        print "No finite data points to use"

