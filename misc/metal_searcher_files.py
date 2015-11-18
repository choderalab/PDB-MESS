from __future__ import print_function
import mdtraj as md
import glob
import multiprocessing as mp

pdbpath = '/cbio/jclab/share/pdb/*/*.ent.gz'
ppn = 32
metals_list = ['TH',
 'BA',
 'YB',
 'EU',
 'FE',
 'DY',
 'V',
 'GA',
 'HG',
 'PR',
 'NI',
 'PT',
 'LA',
 'NA',
 'LI',
 'PB',
 'RE',
 'TL',
 'LU',
 'RU',
 'RB',
 'TE',
 'TB',
 'K',
 'ZN',
 'CO',
 'PD',
 'AG',
 'CA',
 'IR',
 'AM',
 'AL',
 'CE',
 'CD',
 'GD',
 'AU',
 'W',
 'IN',
 'CS',
 'CR',
 'CU',
 'MG',
 'SR',
 'MO',
 'MN',
 'SM',
 'SB',
 'OS',
 'HO']

def metal_searcher(file):
    
    metal_searcher_results = {}
    
    metal_searcher_results['ok_file_count'] = 0
    metal_searcher_results['error_file_count'] = 0
    metal_searcher_results['metals_in_file'] = {}
    metal_searcher_results['any_metal_in_file'] = 0
    
    for i in metals_list:
        metal_searcher_results['metals_in_file'][i] = 0
    
    print(file)
    
    try:
        traj = md.load_pdb(file)
        metal_searcher_results['ok_file_count'] = 1
    except:
        metal_searcher_results['error_file_count'] = 1
        return metal_searcher_results
        
    topo = traj.topology
    
    one_atom_residue_atoms = [atom.name for atom in topo.atoms if atom.name == atom.residue.name and atom.name in metals_list]
    
    if one_atom_residue_atoms:
        metal_searcher_results['any_metal_in_file'] = 1
        for atom_name in one_atom_residue_atoms:
            metal_searcher_results['metals_in_file'][atom_name] += 1
            
    print(metal_searcher_results)        
    return metal_searcher_results

def database_scanner(pdbpath):
    
    database_scanner_results = {}
    
    database_scanner_results['ok_file_count'] = 0
    database_scanner_results['error_file_count'] = 0
    database_scanner_results['metals_in_file'] = {}
    database_scanner_results['any_metal_in_file'] = 0
    
    for i in metals_list:
        database_scanner_results['metals_in_file'][i] = 0
    
    for metal_searcher_results in pool.map(metal_searcher, glob.iglob(pdbpath)):
        database_scanner_results['ok_file_count'] += metal_searcher_results['ok_file_count']
        database_scanner_results['error_file_count'] += metal_searcher_results['error_file_count']
        database_scanner_results['any_metal_in_file'] += metal_searcher_results['any_metal_in_file']
        
        for i in database_scanner_results['metals_in_file']:
            database_scanner_results['metals_in_file'][i] += metal_searcher_results['metals_in_file'][i]
            
    return database_scanner_results        
            
if __name__ == '__main__':
    pool = mp.Pool(processes = ppn)
    database_scanner_results = database_scanner(pdbpath)   
    
with open('metal_searcher_results.txt', 'w') as f:
    f.write(str(database_scanner_results))
    
                 
        
                
