#!/usr/bin/env python
import sys
import copy

mydir = sys.argv[1]
parts = sys.argv[2]
maxx = float(sys.argv[3])/2.
maxy = float(sys.argv[4])/2.

inname  = sys.argv[1]+'/txt_stack/'+parts+'.txt'
outname = sys.argv[1]+'/txt_stack/'+parts+'_ext.txt'
nparts = int(open(inname).readline().rstrip().split('=')[-1])
dat = open(inname, 'r')
out = open(outname, 'w')

prv_traj = [[] for i in range(nparts)]
nbndx = [0 for i in range(nparts)]
nbndy = [0 for i in range(nparts)]

out.write(dat.readline()) #header
fil_index = -1 
for i in range(nparts):
    
    line        = dat.readline()
    coords      = map(float, line.rstrip().split('\t'))
    prv_traj[i] = coords 
    
    coords_adj  = copy.copy(coords)
    
    if coords[3] == fil_index:
        dx = coords[0] - prv_traj[i-1][0]
        dy = coords[1] - prv_traj[i-1][1]
        
        nbndx[i] = nbndx[i-1]
        nbndy[i] = nbndy[i-1]

        if   dx >  maxx: nbndx[i] -= 1
        elif dx < -maxx: nbndx[i] += 1 
        if   dy >  maxy: nbndy[i] -= 1
        elif dy < -maxy: nbndy[i] += 1
        
        coords_adj[0] += 2*nbndx[i]*maxx 
        coords_adj[1] += 2*nbndy[i]*maxy

    out.write('\t'.join(map(str, coords_adj))+'\n')
    fil_index = coords[3]

count = 0

for line in dat:
    
    if count % (nparts+1) == 0:
        out.write(line)
    
    else: 
        n = count % (nparts+1) - 1
        coords = map(float, line.rstrip().split('\t'))

        dx = coords[0] - prv_traj[n][0]
        dy = coords[1] - prv_traj[n][1]
        
        if   dx >  maxx: nbndx[n] -= 1
        elif dx < -maxx: nbndx[n] += 1
        if   dy >  maxy: nbndy[n] -= 1
        elif dy < -maxy: nbndy[n] += 1
        
        coords_adj = copy.copy(coords)
        coords_adj[0] += 2*nbndx[n]*maxx 
        coords_adj[1] += 2*nbndy[n]*maxy
        out.write('\t'.join(map(str, coords_adj))+'\n')
        prv_traj[n] = coords
    
    count += 1

dat.close()
out.close()
