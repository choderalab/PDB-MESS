#!/cbio/jclab/home/rafal.wiewiora/anaconda/bin/python
from __future__ import print_function
import glob
import numpy as np
import mdtraj as md

pdbpath = '../zn/*.pdb'
finaldict = {} 
globfinaldict = {}
progress = 0

for i in glob.iglob(pdbpath):
    try:
        traj = md.load(i)
    except:
        continue

    topo = traj.topology
    bondsindex = topo.to_dataframe()[1]
    if not bondsindex.size:
        continue
    znindex = np.empty((0,2), dtype='int64')
    residues = np.empty((0,), dtype='int64')

    j = 0
    while j < len(bondsindex):
        if topo.atom(bondsindex[j][0]).name == "ZN":
            #znindex = np.append(znindex, [bondsindex[j]], axis=0)
            residues = np.append(residues, [bondsindex[j][1]])
        elif topo.atom(bondsindex[j][1]).name == "ZN":
            #znindex = np.append(znindex, [bondsindex[j]], axis=0)
            residues = np.append(residues, [bondsindex[j][0]])
        j += 1

    if not residues.size:
        continue

    for k in np.nditer(residues):
        if topo.atom(k).residue.name in finaldict:
            finaldict[topo.atom(k).residue.name] += 1
        else:
            finaldict[topo.atom(k).residue.name] = 1

    for key in finaldict:
        if key in globfinaldict:
            globfinaldict[key] += finaldict[key]
        else:
            globfinaldict[key] = finaldict[key]

    print("%d %s" % (progress, i))
    progress += 1

print(globfinaldict)
