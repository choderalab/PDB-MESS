from __future__ import print_function
from __future__ import division
import mdtraj as md
import shutil
import os
import glob

os.chdir("./4")
pdbs = glob.glob("*.pdb")
os.chdir('..')

i = 0
alllistcat = []
C4H0 = []
C3H1 = []
C2H2 = []
C1H3 = []
C0H4 = []
NA = []

while i < len(pdbs):
    path = "./4//%s" % pdbs[i]

    traj = md.load(path)
    topology = traj.topology

    bonds = [bond for bond in topology.bonds]
    bondsstr = [(str(x), str(y)) for x,y in bonds]
    bondsZN = [(x,y) for x,y in bondsstr if "ZN" in x or "ZN" in y]

    listofzincs1 = [x for x,y in bondsZN if "ZN" in x]
    listofzincs2 = [y for x,y in bondsZN if "ZN" in y]
    listofzincs = listofzincs1 + listofzincs2
    listofzincs = list(set(listofzincs))

    dictzincligand = {}
    dictcat = {}
    cyscounter = 0
    hiscounter = 0
    dictcluster = {}
    C4H0t = []
    C3H1t = []
    C2H2t = []
    C1H3t = []
    C0H4t = []
    NAt = []

# build the zinc - ligand dictionary

    for a in listofzincs:
        listofligs1 = [x for x,y in bondsZN if y == a]
        listofligs2 = [y for x,y in bondsZN if x == a]
        listofligs = listofligs1 + listofligs2
        listofligs = list(set(listofligs))
        dictzincligand[a] = listofligs
        # use zing - ligand dictionary to categorise
        for b in listofligs:
            if "CYS" in b:
                cyscounter += 1
            elif "HIS" in b:
                hiscounter += 1
        if cyscounter == 4 and hiscounter == 0:
            dictcat[a] = "C4H0"
            C4H0.append(a)
            C4H0t.append(a)
        elif cyscounter == 3 and hiscounter == 1:
            dictcat[a] = "C3H1"
            C3H1.append(a)
            C3H1t.append(a)
        elif cyscounter == 2 and hiscounter == 2:
            dictcat[a] = "C2H2"
            C2H2.append(a)
            C2H2t.append(a)
        elif cyscounter == 1 and hiscounter == 3:
            dictcat[a] = "C1H3"
            C1H3.append(a)
            C1H3t.append(a)
        elif cyscounter == 0 and hiscounter == 4:
            dictcat[a] = "C0H4"
            C0H4.append(a)
            C0H4t.append(a)
        else:
            dictcat[a] = NA
            NA.append(a)
            NAt.append(a)
        cyscounter = 0
        hiscounter = 0

# dictionary for clustering
    dictcluster['C4H0'] = len(C4H0t)
    dictcluster['C3H1'] = len(C3H1t)
    dictcluster['C2H2'] = len(C2H2t)

    #print("%s: %s" % (pdbs[i], dictcat))

    print("%s: %s" % (pdbs[i], dictcluster))
    i += 1

totalzincs = len(C4H0) + len(C3H1) + len(C2H2) + len(C1H3) + len(C0H4) + len(NA)
print("Stats:")
print("%s: %d" % ("C4H0", len(C4H0)))
print("%s: %d" % ("C3H1", len(C3H1)))
print("%s: %d" % ("C2H2", len(C2H2)))
print("%s: %d" % ("C1H3", len(C1H3)))
print("%s: %d" % ("C0H3", len(C0H4)))
print("%s: %d" % ("NA", len(NA)))
print("%s: %d" % ("Total", totalzincs))

