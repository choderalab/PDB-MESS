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
    metal_ligand_dict_by_cutoff = {}
    metal_ligand_dict_by_CONECT = {}
    ligand_numbers = []
    CONECT_atoms = []
    CUTOFF_atoms = []
    
    dictionary_of_process_counts = {}
    
    # initialize the results dictionary (dictionary_of_process_counts)        
    dictionary_of_process_counts['ok_file_count'] = ok_file_count
    dictionary_of_process_counts['error_file_count'] = error_file_count
    dictionary_of_process_counts['files_with_metal_count'] = files_with_metal_count
    dictionary_of_process_counts['metal_count'] = metal_count
    dictionary_of_process_counts['ligand_numbers'] = ligand_numbers
    dictionary_of_process_counts['CONECT_atoms'] = CONECT_atoms
    dictionary_of_process_counts['CUTOFF_atoms'] = CUTOFF_atoms
    
    
    # load file
    try:
        traj = md.load_pdb(file)
        ok_file_count += 1
        dictionary_of_process_counts['ok_file_count'] = ok_file_count
    except:
        error_file_count += 1
        dictionary_of_process_counts['error_file_count'] = error_file_count
        return dictionary_of_process_counts

    
    # Select atoms and pairs
    topo = traj.topology
    metal_select_name = 'name %s and resname %s' % (metal_name, metal_name)
    metal_atoms = topo.select(metal_select_name)
    metal_all_pairs = topo.select_pairs(metal_select_name, 'symbol O or symbol N or symbol S or symbol Cl')
    
    # No metal - skip, metal - add to count 
    if not metal_atoms.size:
        return dictionary_of_process_counts
    
    files_with_metal_count += 1
    metal_count += len(metal_atoms)
    
    # Prep dictionary with metal indices as keys
    for i in metal_atoms:
        metal_ligand_dict_by_cutoff[i] = []
        metal_ligand_dict_by_CONECT[i] = []
        
    # Compute distances
    try:
        metal_all_distances = md.compute_distances(traj, metal_all_pairs, periodic = False)
    except:
        error_file_count += 1
        dictionary_of_process_counts['error_file_count'] = error_file_count
        return dictionary_of_process_counts   
    
    # Populate the metal_ligand_dict_by_cutoff with atoms closer than cutoff
    for i in range(len(metal_all_distances[0])):
        if metal_all_distances[0, i] < cutoff:
            if metal_all_pairs[i, 0] in metal_atoms:
                metal_ligand_dict_by_cutoff[metal_all_pairs[i, 0]].append(metal_all_pairs[i, 1])
            elif metal_all_pairs[i, 1] in metal_atoms:
                metal_ligand_dict_by_cutoff[metal_all_pairs[i, 1]].append(metal_all_pairs[i, 0])
            
                
    # Look at the CONECT record through topology.bonds
    for i in metal_atoms:
        for bond in topo.bonds:
            if bond[0].index == i:
                metal_ligand_dict_by_CONECT[i].append(bond[1].index)
            elif bond[1].index == i:
                metal_ligand_dict_by_CONECT[i].append(bond[0].index)
                
    # Compare the performance of CONECT and CUTOFF
    for i in metal_atoms:
        ligand_numbers.append((file, i, len(metal_ligand_dict_by_CONECT[i]), len(metal_ligand_dict_by_cutoff[i])))
        
        for j in metal_ligand_dict_by_CONECT[i]:
            CONECT_atoms.append((file, i, str(topo.atom(j))))
        
        for j in metal_ligand_dict_by_cutoff[i]:
            CUTOFF_atoms.append((file, i, str(topo.atom(j))))
        

    # update the results dictionary (ok_file_count and error_file_count have already been done)

    dictionary_of_process_counts['ok_file_count'] = ok_file_count
    dictionary_of_process_counts['error_file_count'] = error_file_count
    dictionary_of_process_counts['files_with_metal_count'] = files_with_metal_count
    dictionary_of_process_counts['metal_count'] = metal_count
    dictionary_of_process_counts['ligand_numbers'] = ligand_numbers
    dictionary_of_process_counts['CONECT_atoms'] = CONECT_atoms
    dictionary_of_process_counts['CUTOFF_atoms'] = CUTOFF_atoms
                
                                                          
    print(dictionary_of_process_counts)
    return dictionary_of_process_counts                                                      
                     
def ligand_scanner_all_database(pdblist):
    
    dictionary_of_database_results = {}
    dictionary_of_database_results['ok_file_count'] = 0
    dictionary_of_database_results['error_file_count'] = 0
    dictionary_of_database_results['files_with_metal_count'] = 0
    dictionary_of_database_results['metal_count'] = 0
    dictionary_of_database_results['ligand_numbers'] = []
    dictionary_of_database_results['CONECT_atoms'] = []
    dictionary_of_database_results['CUTOFF_atoms'] = []
    
    
    for dictionary_of_process_counts in pool.map(ligand_scanner, pdblist):
        for i in dictionary_of_process_counts:
            dictionary_of_database_results[i] += dictionary_of_process_counts[i]
            
    return dictionary_of_database_results
    
# Multiprocess set-up
if __name__ == '__main__':
    
    pool = mp.Pool(processes = ppn)
    dictionary_of_database_results = ligand_scanner_all_database(pdblist)
    

# Write results
f = open('ligand_scanner_results_no_cross_talk.txt', 'w')
f.write(str(dictionary_of_database_results))
f.close()
