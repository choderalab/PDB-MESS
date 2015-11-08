#!/cbio/jclab/home/rafal.wiewiora/anaconda/bin/python
from __future__ import print_function
import mdtraj as md
import glob
import numpy as np

# declare global
pdbpath = '../zn/*.pdb'
cutoff = 0.3
progress = 1
globC4H0count = 0
globC3H1count = 0
globC2H2count = 0
globC1H3count = 0
globC0H4count = 0
globotherscount = 0

for i in glob.iglob(pdbpath):
# init
    try:
        traj = md.load(i)
    except:
        with open("../zn/errors.txt", "rw+") as f:
            f.write("%s\n" % i)
        continue
    topo = traj.topology
    znhiscyspairs = topo.select_pairs('name ZN', '(resn CYS and symbol S) or (resn HIS and symbol N)')
    if not znhiscyspairs.size:
        continue
    distances = md.compute_distances(traj, znhiscyspairs)

# declare per pdb
    j = 0
    k = 0
    znligandpairs = np.empty((0,2), dtype='int32')
    zns = []
    zncodict = {}
    zncysdict = {}
    znhisdict = {}
    pdbC4H0count = 0
    pdbC3H1count = 0
    pdbC2H2count = 0
    pdbC1H3count = 0
    pdbC0H4count = 0
    pdbotherscount = 0

    while j < len(distances[0]):
        if distances[0][j] < cutoff:
            znligandpairs = np.append(znligandpairs, [znhiscyspairs[j]], axis=0)
        j += 1

    while k < len(znligandpairs):
        if topo.atom(znligandpairs[k][0]).name == "ZN":
            zns.append(znligandpairs[k][0])
        if topo.atom(znligandpairs[k][1]).name == "ZN":
            zns.append(znligandpairs[k][1])
        k += 1

    zns = list(set(zns))

# iterate over Zn's
    for x in zns:
        zncodict[x] = 0
        zncysdict[x] = 0
        znhisdict[x] = 0
        l = 0
        while l < len(znligandpairs):
            if znligandpairs[l][0] == x:
                zncodict[x] += 1
                if topo.atom(znligandpairs[l][1]).residue.name == "CYS":
                    zncysdict[x] += 1
                elif topo.atom(znligandpairs[l][1]).residue.name == "HIS":
                    znhisdict[x] += 1
                #extra ligands go here
            if znligandpairs[l][1] == x:
                zncodict[x] += 1
                if topo.atom(znligandpairs[l][0]).residue.name == "CYS":
                    zncysdict[x] += 1
                elif topo.atom(znligandpairs[l][0]).residue.name == "HIS":
                    znhisdict[x] += 1
                #extra ligands ALSO go here
            l += 1
        # update pdb level category counts and database level counts
        if zncodict[x] == 4 and zncysdict[x] == 4 and znhisdict[x] == 0:
            pdbC4H0count += 1
            globC4H0count += 1
        elif zncodict[x] == 4 and zncysdict[x] == 3 and znhisdict[x] == 1:
            pdbC3H1count += 1
            globC3H1count += 1
        elif zncodict[x] == 4 and zncysdict[x] == 2 and znhisdict[x] == 2:
            pdbC2H2count += 1
            globC2H2count += 1
        elif zncodict[x] == 4 and zncysdict[x] == 3 and znhisdict[x] == 1:
            pdbC3H1count += 1
            globC3H1count += 1
        elif zncodict[x] == 4 and zncysdict[x] == 0 and znhisdict[x] == 4:
            pdbC0H4count += 1
            globC0H4count += 1
        else:
            pdbotherscount += 1
            globotherscount += 1

    # print pdb info
    print("%d %s: C4H0 %d, C3H1 %d, C2H2 %d, C1H3 %d, C0H4 %d, others %d" % (progress, i, pdbC4H0count, pdbC3H1count, pdbC2H2count, pdbC1H3count, pdbC0H4count, pdbotherscount))
    progress += 1

# print database info
print("ALL: C4H0 %d, C3H1 %d, C2H2 %d, C1H3 %d, C0H4 %d, others %d" % (globC4H0count, globC3H1count, globC2H2count, globC1H3count, globC0H4count, globotherscount))
