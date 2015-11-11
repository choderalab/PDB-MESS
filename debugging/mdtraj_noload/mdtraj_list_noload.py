from __future__ import print_function
import multiprocessing as mp
import glob
import mdtraj as md

pdbpath = '/cbio/jclab/share/pdb/*/*.ent.gz'
ppn = 32
metal_name = 'ZN'

def loader(file):

    print(file)
    try:
        traj = md.load_pdb(file)
    except:    
        return file
    
def database_analyzer(pdbpath):
    
    with open('mdtraj_list_noload_results', 'w') as f:
        for file in pool.map(loader, glob.iglob(pdbpath)):
            f.write(str(file))
            f.write('\n')
        
    
if __name__ == '__main__':
    
    pool = mp.Pool(processes = ppn)        
    database_analyzer(pdbpath)
