from __future__ import print_function
import mdtraj as md
import glob
import multiprocessing as mp

pdbpath = '/cbio/jclab/share/pdb/*/*.ent.gz'
ppn = 32

def metal_searcher(file):
    
    metal_searcher_results = {}
    metal_searcher_results['one_atom_residue_names'] = []
    metal_searcher_results['ok_file_count'] = 0
    metal_searcher_results['error_file_count'] = 0
    
    print(file)
    
    try:
        traj = md.load_pdb(file)
        metal_searcher_results['ok_file_count'] = 1
    except:
        metal_searcher_results['error_file_count'] = 1
        return metal_searcher_results
        
    topo = traj.topology
    
    one_atom_residue_atoms = [atom for atom in topo.atoms if atom.name == atom.residue.name]
    
    if one_atom_residue_atoms:
        for atom in one_atom_residue_atoms:
            metal_searcher_results['one_atom_residue_names'].append(str(atom.name))
            
    print(metal_searcher_results)        
    return metal_searcher_results

def database_scanner(pdbpath):
    
    database_scanner_results = {}
    database_scanner_results['one_atom_residue_names'] = []
    database_scanner_results['ok_file_count'] = 0
    database_scanner_results['error_file_count'] = 0
    
    for metal_searcher_results in pool.map(metal_searcher, glob.iglob(pdbpath)):
        for i in database_scanner_results:
            database_scanner_results[i] += metal_searcher_results[i]
            
    return database_scanner_results        
            
if __name__ == '__main__':
    pool = mp.Pool(processes = ppn)
    database_scanner_results = database_scanner(pdbpath)   
    
with open('metal_searcher_results.txt', 'w') as f:
    f.write(str(database_scanner_results))
    
                 
        
                
