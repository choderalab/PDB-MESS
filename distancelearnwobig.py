import mdtraj as md
import shutil
import os
import glob
import numpy as np

os.chdir("../pdb")
pdbs = glob.glob("*.pdb")
i = 0
globdistances = np.empty((1,0), dtype='float32')
bigdistances = np.empty((1,0), dtype='float32')
globwobigdistances = np.empty((1,0), dtype='float32')
cutoff = 0.4
k = 0

try:
    while i < len(pdbs):
        traj = md.load(pdbs[i])
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
        print("%d %s %s" % (i ,pdbs[i]), distances)
        i += 1

except:
    shutil.move(pdbs[i], 'errordistancelearn')

# exclude physically irrelevant results, cutoff declared above

while k < len(globdistances):
     if globdistances[k] >= cutoff:
         bigdistances = np.append(bigdistances, [globdistances[k]])
     else:
         globwobigdistances = np.append(globwobigdistances, [globdistances[k]])
     k += 1


print("Mean: %s" % np.mean(globdistances))
print("Median: %s" % np.median(globdistances))
print("Meanw/obig: %s" % np.mean(globwobigdistances))
print("Medianw/obig: %s" % np.median(globwobigdistances))
print("Max: %s" % np.max(globdistances))
print("Min: %s" % np.min(globdistances))
print("Big distances: %s" % bigdistances)
