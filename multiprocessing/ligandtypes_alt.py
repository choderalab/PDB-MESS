from __future__ import print_function
import glob
import numpy as np
import mdtraj as md
import multiprocessing as mp

pdbpath = '*.pdb'

def ligandtypes(i):
    try:
        traj = md.load(i)
    except:
        return None

    topo = traj.topology
    bondsindex = topo.to_dataframe()[1]
    if not bondsindex.size:
        return None
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
        return None

    for k in np.nditer(residues):
        if topo.atom(k).residue.name in resultdict:
            resultdict[topo.atom(k).residue.name] += 1
        else:
            resultdict[topo.atom(k).residue.name] = 1


# make mp workers, execute jobs
if __name__ == '__main__':
    manager = mp.Manager()
    resultdict = manager.dict()
    pool = mp.Pool(processes=16)
    pool.map(ligandtypes, glob.iglob(pdbpath))

# write resultdict to file
f = open('resultligandtypes.txt', 'w')
f.write(str(resultdict))
f.close()
