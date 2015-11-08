from __future__ import print_function
import mdtraj as md
import glob
import numpy as np
import multiprocessing as mp
import json
# declare global
pdbpath = '../zn/*.pdb'
cutoff = 0.3
#progress = 1
#globC4H0count = 0
#globC3H1count = 0
#globC2H2count = 0
#globC1H3count = 0
#globC0H4count = 0
#globotherscount = 0

def ligandcat(i):
# init
    try:
        traj = md.load(i)
    except:
        with open('../zn/errors.txt', 'w') as f:
            f.write("\n%s\n" % i)
        return None

    topo = traj.topology
    znhiscyspairs = topo.select_pairs('name ZN', '(resn CYS and symbol S) or (resn HIS and symbol N)')
    if not znhiscyspairs.size:
        return None
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
            resultdict['C4H0'] += 1
            #globC4H0count += 1
        elif zncodict[x] == 4 and zncysdict[x] == 3 and znhisdict[x] == 1:
            pdbC3H1count += 1
            resultdict['C3H1'] += 1
            #globC3H1count += 1
        elif zncodict[x] == 4 and zncysdict[x] == 2 and znhisdict[x] == 2:
            pdbC2H2count += 1
            resultdict['C2H2'] += 1
            #globC2H2count += 1
        elif zncodict[x] == 4 and zncysdict[x] == 1 and znhisdict[x] == 3:
            pdbC3H1count += 1
            resultdict['C1H3'] += 1
            #globC3H1count += 1
        elif zncodict[x] == 4 and zncysdict[x] == 0 and znhisdict[x] == 4:
            pdbC0H4count += 1
            resultdict['C0H4'] += 1
            #globC0H4count += 1
        else:
            pdbotherscount += 1
            resultdict['Others'] += 1
            otherstuple = (zncodict[x], zncysdict[x], znhisdict[x])
            otherslist.append(otherstuple)

            #globotherscount += 1

# make mp workers, run jobs
if __name__ == '__main__':
    manager = mp.Manager()
    resultdict = manager.dict()
    otherslist = manager.list()
    resultdict['C4H0'] = 0
    resultdict['C3H1'] = 0
    resultdict['C2H2'] = 0
    resultdict['C1H3'] = 0
    resultdict['C0H4'] = 0
    resultdict['Others'] = 0

    pool = mp.Pool(processes=16)
    pool.map(ligandcat, glob.iglob(pdbpath))


# write results to file
#f = open('resultsligandcat.txt', 'w')
#f.write(str(resultdict) + '\n\n\n' + str(othersdictcys) + '\n\n\n' + str(othersdicthis))
#f.close()

# dump dict and list into json
f = open('ligandcatdict.json', 'w')
g = open('ligandcatotherslist.json', 'w')
json.dump(dict(resultdict), f)
json.dump(list(otherslist), g)
f.close()
g.close()
