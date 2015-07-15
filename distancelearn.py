from __future__ import print_function
import mdtraj as md
import shutil
import os
import glob
import numpy as np

os.chdir("./Xray")
pdbs = glob.glob("*.pdb")
os.chdir('..')

i = 0
globdistances = np.empty((1,0), dtype='float32')

while i < len(pdbs):
    path = "./Xray/%s" % pdbs[i]
    traj = md.load(path)
    topology = traj.topology

    allatoms = [atom for atom in topology.atoms]
    allatoms = [str(x) for x in allatoms]
    bondsindex = topology.to_dataframe()[1]
    znindex = np.empty((0,2), dtype='int64')

    j = 0
    while j < len(bondsindex):
        if "ZN" in allatoms[bondsindex[j][0]] or "ZN" in allatoms[bondsindex[j][1]]:
            znindex = np.append(znindex, [bondsindex[j]], axis=0)
        j += 1

    distances = md.compute_distances(traj, znindex)
    globdistances = np.append(globdistances, distances)
    print("%d %s" % (i,pdbs[i]))
    i += 1

print("Mean: %s" % np.mean(globdistances))
print("Median: %s" % np.median(globdistances))
print("Max: %s" % np.max(globdistances))
print("Min: %s" % np.min(globdistances))
