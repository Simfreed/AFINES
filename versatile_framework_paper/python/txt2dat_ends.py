import sys
dr = sys.argv[1]
parts = sys.argv[2]
skip = 1 if len(sys.argv) < 4 else int(sys.argv[3])
tinit = 0 if len(sys.argv) < 5 else int(sys.argv[4])
tfinal = 100000 if len(sys.argv) < 6 else int(sys.argv[5])
nmon = 21 if len(sys.argv) < 7 else int(sys.argv[6])


inname  = dr+'/txt_stack/'+parts+'.txt'
outname = dr+'/txt_stack/'+parts+'.dat'
if len(sys.argv) > 4:
    outname = outname.replace(parts, parts+'{0}-{1}'.format(tinit, tfinal))

dat = open(inname, 'r')
out = open(outname, 'w')

header = dat.readline()
nparts = int(header.rstrip().split('=')[-1])
fil_end_idxs = set(range(1, nparts+1, nmon) + range(nmon, nparts+1, nmon))

tcount = 0
pcount = 0
for line in dat:
    if (line.startswith('t')):
        tcount += 1
        if tcount > tfinal: 
            break
        elif tcount > tinit and tcount % skip == 0:
            out.write('\n')
            pcount = 0
    elif tcount >= tinit and tcount % skip == 0: 
        pcount += 1
        if pcount == nparts:
            out.write('{'+line.rstrip().replace('\t',', ')+'}')
        elif pcount in fil_end_idxs:
            out.write('{'+line.rstrip().replace('\t',', ')+'}\t')

dat.close()
out.close()
