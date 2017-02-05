#!/usr/bin/env python
import sys
import argparse
def write_fibs(ofile, fibdict):
    fibstrs = ['\t'.join(map(lambda x:'{{{0},{1}}}'.format(x[0],x[1]), fibdict[i])) for i in range(len(fibdict))]
    ofile.write('\t'.join(fibstrs)+'\n')

parser = argparse.ArgumentParser()

parser.add_argument("directory", type=str, help="directory of objects.cmo file")
parser.add_argument("--objects", type=str, nargs='+', help="objects to convert", default=['fibers'])

args = parser.parse_args()

infile=open('{0}/objects.cmo'.format(args.directory),'r')
outfile=open('{0}/objects.dat'.format(args.directory), 'w')

read = False
fibs = {}
t = 0
maxt = 50
for line in infile:
    if line.startswith('#section'):
        if line.startswith('#section fiber'): 
            read = True
            fibs = {}
        else:
            read = False
            if line.startswith('#section end'):
                write_fibs(outfile,fibs)
    elif line.startswith('f0:'):
        n = int(line.split(' ')[0].split(':')[-1]) - 1
    elif read:
        xy = map(float,line.split())
        if n in fibs:
            fibs[n].append(xy)
        else:
            fibs[n]=[xy]

outfile.close()
