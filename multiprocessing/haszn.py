from __future__ import print_function
import mdtraj as md
import multiprocessing as mp
import glob

pdbpath = '/cbio/jclab/share/pdb/*/*.ent.gz'

def haszn(i):
    try:
        traj = md.load_pdb(i)
    except:
        return None

    topo = traj.topology

    atoms = [str(x) for x in topo.atoms]
    atomszn = [x for x in atoms if "ZN" in x]

    if atomszn:
        znpdbs.append(i)

# multiprocessing
if __name__ == '__main__':
    manager = mp.Manager()
    znpdbs = manager.list()
    pool = mp.Pool(processes = 16)
    pool.map(haszn, glob.iglob(pdbpath))

# get total pdb files no.
alldatabase = glob.glob(pdbpath)

# write to file
f = open('hasznresults.txt', 'w')
f.write(str(len(znpdbs)))
f.write(str(len(alldatabase)))
f.close()
