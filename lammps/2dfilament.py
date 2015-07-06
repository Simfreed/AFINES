import math
import random
import sys
import time

xbds=[-500,500]
ybds=[-500,500]
zbds=[-0.05,0.05]

bondlength = float(sys.argv[1]) if len(sys.argv) > 0 else 1
nrod = int(200 / bondlength)

seed = sys.argv[2] if len(sys.argv) > 1 else time.time()
random.seed(seed)

rod_pos=[]

pos = [-150,0,0]
phi = 0 
small_angle = 0.1 #math.pi/12
nbondtypes = 3 
#atoms
atoms = []

variance = 0.1 #kT/(bending_stiffness*bondLength*bondLength)
for i in range(nrod):
    
    rod_pos.append(pos)

    phi += random.normalvariate(0, variance)

    atom   = [i+1, 1, 1, pos[0], pos[1], 0]#, pos[0] + rodlength*math.cos(phi), pos[0] + rodlength*math.sin(phi), 0]
    atoms.append(atom)

    pos = (pos[0] + bondlength*math.cos(phi), pos[1]+bondlength*math.sin(phi), 0)

#bonds
bonds = []
for i in range(nrod-1):
    
    bonds.append([i, i%nbondtypes + 1, i+1, i+2])

#angles
angles = []
for i in range(nrod-2):

    angles.append([i+1, 1, i+1, i+2, i+3])

#write to file
f=open('init_cfgs/filament_bl{0:.2f}.txt'.format(bondlength),'w')
f.write('#Actin filament\n')
f.write('\n{0} atoms\n{1} bonds\n{2} angles\n'.format(len(atoms),len(bonds),len(angles)))
f.write('\n1 atom types\n{0} bond types\n1 angle types\n'.format(nbondtypes))

f.write('\n{0} {1} xlo xhi\n{2} {3} ylo yhi\n{4} {5} zlo zhi\n'.format(xbds[0], xbds[1], ybds[0], ybds[1], zbds[0],
    zbds[1]))

f.write('\nMasses\n\n1 1\n')

f.write('\n\nAtoms\n\n')
for a in atoms:
    f.write(' '.join([str(i) for i in a])+'\n')

f.write('\n\nBonds\n\n')
for a in bonds:
    f.write(' '.join([str(i) for i in a])+'\n')

f.write('\n\nAngles\n\n')
for a in angles:
    f.write(' '.join([str(i) for i in a])+'\n')

f.close()
