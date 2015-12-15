from __future__ import print_function
import mdtraj as md
import glob
import multiprocessing as mp

pdblist_file = open('pdblist.txt')
pdblist = [line[:-1] for line in pdblist_file]
cutoff = 0.3
metal_name = 'ZN'
ppn = 32

def ligand_scanner(file):
    
    print(file)
    
    ok_file_count, error_file_count, metal_count, files_with_metal_count = 0, 0, 0, 0
    metal_ligand_dict_by_CONECT = {}
    dictionary_of_process_counts = {}
    ligand_distances = []

    
    # initialize the results dictionary (dictionary_of_process_counts)        
    dictionary_of_process_counts['ok_file_count'] = ok_file_count
    dictionary_of_process_counts['error_file_count'] = error_file_count
    #dictionary_of_process_counts['files_with_metal_count'] = files_with_metal_count
    #dictionary_of_process_counts['metal_count'] = metal_count
    dictionary_of_process_counts['ligand_distances'] = ligand_distances
    
    # load file
    try:
        traj = md.load_pdb(file)
        ok_file_count += 1
        dictionary_of_process_counts['ok_file_count'] = ok_file_count
    except:
        error_file_count += 1
        dictionary_of_process_counts['error_file_count'] = error_file_count
        return dictionary_of_process_counts
        
    topo = traj.topology
    metal_select_name = 'name %s and resname %s' % (metal_name, metal_name)
    metal_atoms = topo.select(metal_select_name)
    
    cys_ligand_atoms = [atom.index for atom in topo.atoms if atom.name == 'SG' and atom.residue.name == 'CYS']
    his_ligand_atoms = [atom.index for atom in topo.atoms if (atom.name == 'ND1' or atom.name == 'NE2') and atom.residue.name == 'HIS']
    glu_ligand_atoms = [atom.index for atom in topo.atoms if (atom.name == 'OE1' or 'OE2') and atom.residue.name == 'GLU']
    asp_ligand_atoms = [atom.index for atom in topo.atoms if (atom.name == 'OD1' or 'OD2') and atom.residue.name == 'ASP']
    
    ligand_atoms = cys_ligand_atoms + his_ligand_atoms + glu_ligand_atoms + asp_ligand_atoms
    
    for i in metal_atoms:
        metal_ligand_dict_by_CONECT[i] = []
        
    # Look at the CONECT record through topology.bonds
    for i in metal_atoms:
        for bond in topo.bonds:
            if bond[0].index == i:
                metal_ligand_dict_by_CONECT[i].append(bond[1].index)
            elif bond[1].index == i:
                metal_ligand_dict_by_CONECT[i].append(bond[0].index)    
                
                
    # Analyze distances in CONECT
    for i in metal_ligand_dict_by_CONECT:
        ligands_for_distance_comp = [j for j in metal_ligand_dict_by_CONECT[i] if j in ligand_atoms]

        metal_all_pairs = topo.select_pairs([i], ligands_for_distance_comp)
    
        try:
            metal_all_distances = md.compute_distances(traj, metal_all_pairs, periodic = False) 
            for k in metal_all_distances:
                ligand_distances.append(k) 
        except:
            pass         
            
          
            
             
    # Update dictionary of resultsdictionary_of_process_counts['ok_file_count'] = ok_file_count
    dictionary_of_process_counts['ok_file_count'] = ok_file_count
    dictionary_of_process_counts['error_file_count'] = error_file_count
    #dictionary_of_process_counts['files_with_metal_count'] = files_with_metal_count
    #dictionary_of_process_counts['metal_count'] = metal_count
    dictionary_of_process_counts['ligand_distances'] = ligand_distances        
    
    return dictionary_of_process_counts
    
def ligand_scanner_all_database(pdblist):
    
    dictionary_of_database_results = {}
    dictionary_of_database_results['ok_file_count'] = 0
    dictionary_of_database_results['error_file_count'] = 0
    #dictionary_of_database_results['files_with_metal_count'] = 0
    #dictionary_of_database_results['metal_count'] = 0
    dictionary_of_database_results['ligand_distances'] = []
    
    for dictionary_of_process_counts in pool.map(ligand_scanner, pdblist):
        for i in dictionary_of_process_counts:
            dictionary_of_database_results[i] += dictionary_of_process_counts[i]
            
    return dictionary_of_database_results
    
# Multiprocess set-up
if __name__ == '__main__':
    
    pool = mp.Pool(processes = ppn)
    dictionary_of_database_results = ligand_scanner_all_database(pdblist)
    

# Write results
f = open('24Nov_dummy_atom_scheme_conect_distances.txt', 'w')
f.write(str(dictionary_of_database_results))
f.close()        
                
                    
                        
