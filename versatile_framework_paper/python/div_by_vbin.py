#!/usr/bin/env python
import numpy as np
from scipy.interpolate import Rbf
from scipy.stats import binned_statistic_2d
import os
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("directory", type=str, help="directory of txt_stack/actins.txt file")
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

fname = args.directory + '/txt_stack/actins.txt'

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

xi = np.arange(-0.5*xsz, 0.5*xsz, args.dr[0])
yi = np.arange(-0.5*ysz, 0.5*ysz, args.dr[1])

divxydir = '{0}/analysis/xyuvd_nbins{1}_lsc{2}'.format(args.directory,args.nbins,args.div_lsc) 
div_bins = []
if not os.path.exists(divxydir):
    os.makedirs(divxydir)

f      = 0.5
edgesx = np.arange(-f*xsz, f*xsz + args.div_lsc, args.div_lsc)
edgesy = np.arange(-f*ysz, f*ysz + args.div_lsc, args.div_lsc)

for t in range(len(mxyt)):

    idxs = np.logical_and( np.isfinite( mxyt[t,:,0] ) , np.isfinite(mxyt[t,:,1]) )
    
    if len(mxyt[t,idxs]) > 0:
        
        rbfu = Rbf(mxyt[t,idxs,0], mxyt[t,idxs,1], muvt[t,idxs,0], function=args.rbfunc, epsilon=args.rbeps)
        rbfv = Rbf(mxyt[t,idxs,0], mxyt[t,idxs,1], muvt[t,idxs,1], function=args.rbfunc, epsilon=args.rbeps)

        # set the grid points (using midpoints because of totdiv calculation)
        xgrid, ygrid = np.meshgrid((xi[1:]+xi[:-1])*0.5, (yi[1:]+yi[:-1])*0.5)

        # evaluate interpolation on grid
        uint = rbfu(xgrid, ygrid)
        vint = rbfv(xgrid, ygrid)

        # calculate the divergence
        # first element in gradient list is gradient of rows (i.e.,  d/dy)
        # second '''''''''''''''''''''''''''''''''''''' cols (i.e.,  d/dx)
        divxy  = np.gradient(uint, args.dr[0])[1] + np.gradient(vint, args.dr[1])[0] * args.dr[0] * args.dr[1];
        xyuvd  = np.stack([xgrid.flatten(), ygrid.flatten(), uint.flatten(), vint.flatten(), divxy.flatten()], axis=1)
        #print 'xgrid.flatten().shape ={0}, ygrid.flatten.shape()={1},divxy.flatten.shape()={2}'.format(xgrid.flatten().shape, ygrid.flatten().shape, divxy.flatten().shape)
        div_bins.append(binned_statistic_2d(xgrid.flatten(), ygrid.flatten(), divxy.flatten(), statistic='sum',
            bins=[edgesx, edgesy])[0].flatten())
         
        np.savetxt('{0}/t{1}.txt'.format(divxydir, t), xyuvd, '%10.5f')
        #print 't={0}\tall divergences={1}'.format(t, div_bins[t])
        print 't = {0}\tmin divergence = {1}'.format(t, np.min(div_bins[t]))
    else:
        print "No finite data points to use"

np.savetxt('{0}/analysis/div_bins_nbins{1}_lsc{2}.txt'.format(args.directory, args.nbins, args.div_lsc), np.array(div_bins), '%10.5f')
#    divtott.append(np.sum(divxy)*dx*dy)
#    xyuvdt.append( np.stack([xgrid, ygrid, uint, vint, divxy], axis=1) )


