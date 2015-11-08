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
    znconectlist = []
    haszn = False
    global znpdbs
    global allinplace
    global listofcoordnos
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
                        znconectlist.append(zn)
                        listofcoordnos.append(len(fields) - 2)
    if haszn == True and len(set(listofzn)) == len(set(znconectlist)):
         with lock:
             allinplace.value += 1
            
            

# multiprocessing
if __name__ == '__main__':
    lock = mp.Lock()
    manager = mp.Manager()
    znpdbs = mp.Value(c_int)
    allinplace = mp.Value(c_int)
    listofcoordnos = manager.list()
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
f.write("The coordination list:\n")
f.write(str(listofcoordnos))
f.close()
