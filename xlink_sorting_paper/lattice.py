import numpy as np
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("directory", type=str, help="output directory")
parser.add_argument("--dL", type=float, help="difference in xlink length (um)", default = 0.1)
parser.add_argument("--nSites", type=int, help="number of starting lattice sites", default=405)
parser.add_argument("--nSitesMax", type=int, help="max number of lattice sites", default=405)
parser.add_argument("--muA", type=float, help="chemical potential of spacer A (pNum)", default=0.004)
parser.add_argument("--muB", type=float, help="chemical potential of spacer B (pNum)", default=0.004)
parser.add_argument("--tf", type=int, help="number of time steps",default=pow(10,7))
parser.add_argument("--kT", type=float, help="temperature*(boltzmann's constant) (pNum)", default=0.004)
parser.add_argument("--Lp", type=float, help="filament persistence length (um)", default=17)
parser.add_argument("--lSite", type=float, help="length of lattice site (um)", default=0.037)
parser.add_argument("--seed", type=int, help="seed for random number generator", default=None)
parser.add_argument("--checkEng", type=bool, help="flag for checking energy", default=False)
parser.add_argument("--growthPeriod", type=int, help="add a site every nSites[t]*growthPeriod steps", default=-1)

parser.add_argument("--muAnn", type=float, help="cooperativity for nearest neighbor interaction for spacer A (pNum)", default=0)
parser.add_argument("--muBnn", type=float, help="cooperativity for nearest neighbor interaction for spacer B (pNum)", default=0)

parser.add_argument("--dumpPeriod", type=int, help="number of time steps between state dumps", default=1000)
parser.add_argument("--dumpAfter", type=int, help="minimum iteration on which to dump", default=-1)

args = parser.parse_args()
with open('{0}/args.txt'.format(args.directory),'w') as argfile: 
    for arg in vars(args):
        argfile.write('{0}={1}\n'.format(arg, getattr(args, arg)))

tol = 0.00001

np.seterr(invalid='print')
np.random.seed(args.seed)
kappaB = args.Lp*args.kT
gapLengths = np.arange(1, args.nSitesMax)*args.lSite
gapEngs = kappaB / gapLengths*np.power(np.arcsin(args.dL/gapLengths),2)
#print(gapEngs)

infty=np.nanmax(gapEngs)*10000
gapEngs[np.where(gapLengths<=args.dL)]=infty
#print(gapEngs)

# 0 in -1th index is for non-existing lattice site (i.e., hasn't grown to that length yet)
mus = [0, args.muA, args.muB, 0] 
musnn = [0, args.muAnn, args.muBnn, 0] 

def get_eng(lattice):
    eng = 0
    gapStarter = 0
    nGap = 0 
    prev = -1
    for i in lattice:
        
        if i == -1: # end of lattice
            continue

        eng -= mus[i]
        
        # constant cooperativity
        if i == prev:
            eng -= musnn[i]
        prev = i

        # gap calculation
        if i == 0:
            nGap += 1
        else:
            if gapStarter != i and gapStarter != 0:
                eng += gapEngs[nGap]
            nGap = 0
            gapStarter = i
        

    return eng

def metropolis_prob(dEng):
    return min(1, np.exp(-dEng/args.kT))

def two_gap_eng(leftState, middleState, rightState, leftGap, rightGap):

    if middleState==0:
        return gap_eng(leftState, rightState, leftGap+rightGap+1)
    else:
        return gap_eng(leftState, middleState, leftGap) + gap_eng(middleState, rightState, rightGap)

def gap_eng(leftState, rightState, nGap):
    if leftState == rightState or leftState == -1 or rightState == -1: # -1 means end of chain
        return 0
    else:
        return gapEngs[nGap]

def get_dEng(lattice, site, newState):
    
    oldState = lattice[site]
    if oldState == newState:
        return 0
    
    else:
        nRightGap = 0
        rightBound = -1
        nLeftGap = 0
        leftBound = -1
        
        muCoop = 0
        
        for i in range(site+1, args.nSitesMax):
            if lattice[i] == 0:
                nRightGap += 1
            else:
                rightBound = lattice[i]
                break
        
        for i in range(site-1, -1, -1):
            if lattice[i] == 0:
                nLeftGap += 1
            else:
                leftBound = lattice[i]
                break
        
        rightState  = lattice[site+1] if site + 1 < args.nSites else -1
        leftState   = lattice[site-1] if site - 1 >= 0 else -1
        muCoop     -= musnn[newState] * ((newState == rightState) + (newState == leftState))
        muCoop     += musnn[oldState] * ((oldState == rightState) + (oldState == leftState))
       
        return  - mus[newState] + mus[oldState]\
                + two_gap_eng(leftBound, newState, rightBound, nLeftGap, nRightGap)\
                - two_gap_eng(leftBound, oldState, rightBound, nLeftGap, nRightGap)\
                + muCoop
       

nSitesT     = args.nSites
lattice     = [0]*nSitesT + [-1]*(args.nSitesMax - nSitesT)
#eng = get_eng(lattice)
eng = 0

latticeOut  = open('{0}/lattice.dat'.format(args.directory), 'w')
engOut      = open('{0}/engs.dat'.format(args.directory), 'w')
gT          = args.growthPeriod if args.growthPeriod > 0 else args.tf + 1
lastGrowT   = -1 
for t in range(args.tf):

    if t > args.dumpAfter and t % args.dumpPeriod == 0:
        latticeOut.write(','.join(map(str,lattice))+'\n')
        engOut.write('{0}\n'.format(eng))

    site = np.random.randint(nSitesT)
    state = np.random.randint(3)
    
    #nlattice = lattice 
    #nlattice[site] = state
    #neng = get_eng(nlattice)
    deng = get_dEng(lattice, site, state)

    if np.random.random() < metropolis_prob(deng):
        lattice[site] = state
        eng += deng
        if args.checkEng:
            cEng = get_eng(lattice)
            if cEng != 0 and np.fabs(eng - cEng)/cEng > tol:
                print("\nERROR: ENERGY = {0} DOESN'T MATCH PREDICTED ENERGY = {1} WITHIN TOLERANCE OF {2}%".format(eng,
                    cEng, tol*100))

    if (t-lastGrowT) % (gT*nSitesT) == 0 and nSitesT < args.nSitesMax:
        nSitesT += 1
        lattice[nSitesT-1]=0
        lastGrowT = t

latticeOut.close()
engOut.close()
