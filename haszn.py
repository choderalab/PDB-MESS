from __future__ import print_function
import mdtraj as md
import shutil

for i in glob.iglob('*.pdb'):
    traj = md.load(i)
    topo = traj.topology

    atoms = [x for x in topo.atoms]
    atomszn = [x for x in atoms if "ZN" in x]

    if atomszn:
        shutil.move(i, 'zn')

    print(i)
