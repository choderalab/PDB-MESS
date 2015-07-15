import mdtraj as md
import shutil
import os
import glob

os.chdir("./pdbs")
pdbs = glob.glob("*.pdb")
os.chdir('..')

i = 0
while i < len(pdbs):
    path = "./pdbs//%s" % pdbs[i]
    pathconnect = "./connect/%s" % pdbs[i]
    patherror = "./pdbserror/%s" % pdbs[i]

    try:
        traj = md.load(path)
        topology = traj.topology

        bonds = [bond for bond in topology.bonds]
        bondsstr = [(str(x), str(y)) for x,y in bonds]
        bondsZN = [(x,y) for x,y in bondsstr if "ZN" in x or "ZN" in y]

        if bondsZN:
            shutil.move(path, pathconnect)
            print(pdbs[i])

    except:
        shutil.move(path, patherror)
        print("%s ERROREROOREROOR" % pdbs[i])

    i += 1

