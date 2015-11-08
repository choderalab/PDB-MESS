from __future__ import print_function
import multiprocessing as mp
import glob
import gzip
from ctypes import c_int

pdbpath = '/cbio/jclab/share/pdb/*/*.ent.gz'
processes = 16

def znconect(i):
    pdbfile = gzip.open(i)
    listofzn = []
    znconectlogicdict = {}
    haszn = False
    #znconectligandnolist = []
    global znpdbs
    global allinplace
    for line in pdbfile:
        fields = line.split()
        if haszn == False:
            if fields[0] == "HET":
                if fields[1] == "ZN":
                    haszn = True
                    with lock:
                        znpdbs.value += 1
        if haszn == True:
            if fields[0] == "HETATM":
                if fields[2] == "ZN":
                    listofzn.append(fields[1])
            elif fields[0] == "CONECT":
                for zn in listofzn:
                    if fields[1] == zn:
                        znconectlogicdict[zn] = True
                        #znconectligandnolist.append(len(line.split()) - 2)
                        #ligandnos.append(len(list.split()) - 2)
    if haszn == True and len(listofzn) == len(znconectlogicdict):
         with lock:
             allinplace.value += 1

# multiprocessing
if __name__ == '__main__':
    #manager = mp.Manager
    lock = mp.Lock()
    znpdbs = mp.Value(c_int)
    allinplace = mp.Value(c_int)
    #ligandnos = manager.list()
    #somemissing = manager.list()
    pool = mp.Pool(processes)
    pool.map(znconect, glob.iglob(pdbpath))

# get total pdb files no.
alldatabase = glob.glob(pdbpath)

# write to file
f = open('znconectresults', 'w')
f.write("All Zn-cont. pdb's:\n")
f.write(str(znpdbs.value))
f.write("\n")
f.write("PDB's with complete CONECT record:\n")
f.write(str(allinplace.value))
f.write("\n")
f.write("All files:\n")
f.write(str(len(alldatabase)))
f.write("\n")
#f.write("Ligand nos:\n")
#f.write(ligandnos)
f.close()
