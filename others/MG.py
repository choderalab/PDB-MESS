from __future__ import print_function
import multiprocessing as mp
import glob
import gzip
from ctypes import c_int
import json

pdbpath = '/cbio/jclab/share/pdb/*/*.ent.gz'
processes = 32

def metalconect(i):
    pdbfile = gzip.open(i)
    listofmetals = []
    metalconectlist = []
    hasmetal = False
    ligandsdict = {}
    global metalpdbs
    global allinplace
    global listofcoordnos
    global highcooordnospaths
    for line in pdbfile:
        fields = line.split()
        if hasmetal == False:
            if fields[0] == "HET":
                if fields[1] == "MG":
                    hasmetal = True
                    with lock:
                        metalpdbs.value += 1
        if hasmetal == True:
            if fields[0] == "HETATM":
                if fields[2] == "MG":
                    listofmetals.append(fields[1])                
    
    listofmetals = list(set(listofmetals))
    for metal in listofmetals:
        ligandsdict[metal] = 0
        
    if hasmetal == True:
        pdbfile.seek(0)
        for line in pdbfile:
            fields = line.split()
            if fields[0] == "CONECT":
                for metal in listofmetals:
                    if fields[1] == metal:
                        metalconectlist.append(metal)
                        ligandsdict[metal] += len(fields) - 2     
        
    
    if hasmetal == True and len(set(listofmetals)) == len(set(metalconectlist)):
        with lock:
             allinplace.value += 1
             
    for metal in ligandsdict:
        if ligandsdict[metal] > 10:
            highcoordnospaths[i] = ligandsdict[metal]
            listofcoordnos.append(ligandsdict[metal])         
        else:
            listofcoordnos.append(ligandsdict[metal])    
    
    print(i)        
    pdbfile.close()
    
# multiprocessing
if __name__ == '__main__':
    lock = mp.Lock()
    manager = mp.Manager()
    metalpdbs = mp.Value(c_int)
    allinplace = mp.Value(c_int)
    listofcoordnos = manager.list()
    highcoordnospaths = manager.dict()
    pool = mp.Pool(processes)
    pool.map(metalconect, glob.iglob(pdbpath))

# get total pdb files no.
alldatabase = glob.glob(pdbpath)

# write to file
f = open('metalresults', 'w')
f.write("All metal-cont. pdb's:\n")
f.write(str(metalpdbs.value))
f.write("\n")
f.write("PDB's with complete CONECT record:\n")
f.write(str(allinplace.value))
f.write("\n")
f.write("All files:\n")
f.write(str(len(alldatabase)))
f.write("\n")
f.write("ligand coords list:\n")
f.write(str(listofcoordnos))
f.write("\n")
f.write("coord nos over 10:\n")
f.write(str(highcoordnospaths))
f.close()
# write listofcoordnos to json 
#json.dump(listofcoordnos, open('listofcoordnos.json', 'w'))
    
