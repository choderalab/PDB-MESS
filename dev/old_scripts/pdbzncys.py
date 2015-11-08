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
    pathZN = "./pdbsZNCys/%s" % pdbs[i]
    patherror = "./pdbserror/%s" % pdbs[i]
    
    try:
        traj = md.load(path)
        topology = traj.topology

        bonds = [bond for bond in topology.bonds]
        bondsstr = [(str(x), str(y)) for x,y in bonds]
        bondsZN = [(x,y) for x,y in bondsstr if ("ZN" in x and "CYS" in y) or ("ZN" in y and "CYS" in x)]

        if bondsZN:
            shutil.move(path, pathZN)

    except:
        shutil.move(path, patherror)

    i += 1

