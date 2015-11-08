from __future__ import print_function
import multiprocessing as mp
import glob
import gzip

pdbpath = '/cbio/jclab/share/pdb/*/*.ent.gz'
processes = 16

def znhet(i):
    pdbfile = gzip.open(i)
    for line in pdbfile:
        fields = line.split()
        if fields[0] == "HET":
            if fields[1] == "ZN":
                znpdbs.append(i)
                break
    pdbfile.close()

# multiprocessing
if __name__ == '__main__':
    manager = mp.Manager()
    znpdbs = manager.list()
    pool = mp.Pool(processes)
    pool.map(znhet, glob.iglob(pdbpath))

# get total pdb files no.
alldatabase = glob.glob(pdbpath)

# write to file
f = open('znhetresults', 'w')
f.write("Zn centers found:\n")
f.write(str(len(znpdbs)))
f.write("\n")
f.write("All files:\n")
f.write(str(len(alldatabase)))
f.close()
