import mdtraj as md
#import glob 
import multiprocessing as mp 

#files = glob.iglob('/Users/rafalpwiewiora/pdb/pdb/*/*.ent.gz')
pdblist_file = open('pdblist.txt')
pdblist = [line[:-1] for line in pdblist_file]

def loader(file):
    try:
        traj = md.load_pdb(file)
    except:
        pass
            
pool = mp.Pool(processes = 8)
pool.map(loader, pdblist)    
    
    
