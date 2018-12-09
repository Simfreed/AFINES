import sys
import numpy as np
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("outfile", type=str, help="output")
parser.add_argument("--fov", type=float, nargs=2, help="field of view", default=[35,20])
parser.add_argument("--lm", type=float, help="rest length of mot", default=0.5)
parser.add_argument("--nm", type=int, help="number of motors", default=10)
parser.add_argument("--seed", type=int, help="random number seed", default=1)

'''
create a list of bins
loop through the bins, finding points in them using np.logical_and
calculate the divergence using rbf interpolation only in that box
'''

args = parser.parse_args()

np.random.seed(args.seed)

x0s = np.random.uniform(-args.fov[0]/2.,args.fov[0]/2.,args.nm)
y0s = np.random.uniform(-args.fov[1]/2.,args.fov[1]/2.,args.nm)
ths = np.random.uniform(-np.pi,np.pi,args.nm)
out = np.vstack([x0s, y0s, args.lm*np.cos(ths), args.lm*np.sin(ths), -1*np.ones(args.nm), -1*np.ones(args.nm), -1*np.ones(args.nm),
    -1*np.ones(args.nm)]).transpose()
np.savetxt(args.outfile, out, fmt="%.5f\t%.5f\t%.5f\t%.5f\t%d\t%d\t%d\t%d")
