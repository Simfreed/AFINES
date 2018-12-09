import sys
import numpy as np
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("outfile", type=str, help="output")
parser.add_argument("--r0", type=float, nargs=2, help="position of head 0 of motor 0", default=[-12.5,-0.25])
parser.add_argument("--dr", type=float, nargs=2, help="distance between mots", default=[0.1,0])
parser.add_argument("--lm", type=float, help="length of mot", default=0.5)
parser.add_argument("--thm", type=float, help="angle of mot", default=np.pi/2)
parser.add_argument("--nm", type=int, help="number of motors", default=10)

'''
create a list of bins
loop through the bins, finding points in them using np.logical_and
calculate the divergence using rbf interpolation only in that box
'''

args = parser.parse_args()

rest_of_pos = [args.lm*np.cos(args.thm), args.lm*np.sin(args.thm),0,1,0,0]
r=args.r0
out = open(args.outfile, 'w')
for i in range(args.nm):
    pos = r+rest_of_pos
    out.write('\t'.join(map(str,pos))+'\n')
    r=[r[0]+args.dr[0], r[1]+args.dr[1]]

out.close()
