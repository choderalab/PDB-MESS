from __future__ import print_function
import mdtraj as md
import shutil
import os
import glob

os.chdir("./Xray")
pdbs = glob.glob("*.pdb")
os.chdir('..')

i = 0
globcoordlist = []
globdictzinc = {}

while i < len(pdbs):
    path = "./Xray//%s" % pdbs[i]

    traj = md.load(path)
    topology = traj.topology

    bonds = [bond for bond in topology.bonds]
    bondsstr = [(str(x), str(y)) for x,y in bonds]
    bondsZN = [(x,y) for x,y in bondsstr if "ZN" in x or "ZN" in y]

    listofzincs1 = [x for x,y in bondsZN if "ZN" in x]
    listofzincs2 = [y for x,y in bondsZN if "ZN" in y]
    listofzincs = listofzincs1 + listofzincs2
    listofzincs = list(set(listofzincs))

    tempdictzinc = {}

    for a in listofzincs:
        listofcys1 = [x for x,y in bondsZN if y == a]
        listofcys2 = [y for x,y in bondsZN if x == a]
        listofcys = listofcys1 + listofcys2
        listofcys = list(set(listofcys))
        globdictzinc[a] = listofcys
        globcoordlist.append(len(listofcys))
        tempdictzinc[a] = len(listofcys)


    print("%s: %s" % (pdbs[i], tempdictzinc))
    i += 1

globcoordlist = list(set(globcoordlist))
print("All coord. numbers: %s" % globcoordlist)    

