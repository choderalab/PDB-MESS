from __future__ import print_function
import mdtraj as md
import glob
import multiprocessing as mp

pdbpath = '/cbio/jclab/share/pdb/*/*.ent.gz'
cutoff = 0.3
metal_name = 'ZN'
ppn = 32


def ligand_scanner(file):
    
    #ok_file_count, error_file_count, metal_count, files_with_metal_count = 0, 0, 0, 0
    #extra_in_cutoff_count_by_atom = 0
    #extra_in_cutoff_count_by_file = 0
    #extra_in_cutoff_file_list = []
    #extra_in_cutoff_residue_list = []
    
    metal_ligand_dict_by_cutoff = {}
    metal_ligand_dict_by_CONECT = {}
    metal_ligand_dict_in_cutoff_not_CONECT = {}
    metal_ligand_dict_CORRECTED_in_cutoff_not_CONECT = {}
    metal_ligand_dict_in_CONECT_not_cutoff = {}
    #extra_in_cutoff_for_file_analyzed = False
    #dictionary_of_process_counts = {}
    
    # initialize the results dictionary (dictionary_of_process_counts)        
    #dictionary_of_process_counts['ok_file_count'] = ok_file_count
    #dictionary_of_process_counts['error_file_count'] = error_file_count
    #dictionary_of_process_counts['files_with_metal_count'] = files_with_metal_count
    #dictionary_of_process_counts['metal_count'] = metal_count
    #dictionary_of_process_counts['extra_in_cutoff_count_by_atom'] = extra_in_cutoff_count_by_atom
    #dictionary_of_process_counts['extra_in_cutoff_count_by_file'] = extra_in_cutoff_count_by_file
    #dictionary_of_process_counts['extra_in_cutoff_file_list'] = extra_in_cutoff_file_list
    #dictionary_of_process_counts['extra_in_cutoff_residue_list'] = extra_in_cutoff_residue_list
    
    numbers_of_extra_cutoff_ligands = []

    # load file
    try:
        traj = md.load_pdb(file)
        #ok_file_count += 1
        #dictionary_of_process_counts['ok_file_count'] = ok_file_count
    except:
        #error_file_count += 1
        #dictionary_of_process_counts['error_file_count'] = error_file_count
        return None

    
    # Select atoms and pairs
    topo = traj.topology
    metal_select_name = 'name %s' % metal_name
    metal_atoms = topo.select(metal_select_name)
    metal_all_pairs = topo.select_pairs(metal_select_name, 'all')
    
    # No metal - skip, metal - add to count 
    if not metal_atoms.size:
        return None
    
    #files_with_metal_count += 1
    #metal_count += len(metal_atoms)
    
    # Prep dictionary with metal indices as keys
    for i in metal_atoms:
        metal_ligand_dict_by_cutoff[i] = []
        metal_ligand_dict_by_CONECT[i] = []
        metal_ligand_dict_in_cutoff_not_CONECT[i] = []
        metal_ligand_dict_in_CONECT_not_cutoff[i] = []
        metal_ligand_dict_CORRECTED_in_cutoff_not_CONECT[i] = []
        
        
    # Compute distances
    metal_all_distances = md.compute_distances(traj, metal_all_pairs)
    
    # Populate the metal_ligand_dict_by_cutoff with atoms closer than cutoff
    for i in range(len(metal_all_distances[0])):
        if metal_all_distances[0, i] < cutoff:
            if metal_all_pairs[i, 0] in metal_atoms:
                metal_ligand_dict_by_cutoff[metal_all_pairs[i, 0]].append(metal_all_pairs[i, 1])
            elif metal_all_pairs[i, 1] in metal_atoms:
                metal_ligand_dict_by_cutoff[metal_all_pairs[i, 1]].append(metal_all_pairs[i, 0])
                
    # Look at the CONECT record through topology.bonds
    for i in metal_ligand_dict_by_CONECT:
        for bond in topo.bonds:
            if bond[0].index == i:
                metal_ligand_dict_by_CONECT[i].append(bond[1].index)
            elif bond[1].index == i:
                metal_ligand_dict_by_CONECT[i].append(bond[0].index)
                
    # Compare the contents of metal_ligand_dict_by_cutoff and metal_ligand_dict_by_CONECT, correct for cutoff recognizing extra atoms for residues already in CONECT
    for i in metal_atoms:
        metal_ligand_dict_in_cutoff_not_CONECT[i] = [x for x in metal_ligand_dict_by_cutoff[i] if x not in metal_ligand_dict_by_CONECT[i]]
        metal_ligand_dict_in_CONECT_not_cutoff[i] = [x for x in metal_ligand_dict_by_CONECT[i] if x not in metal_ligand_dict_by_cutoff[i]]    
    
    for i in metal_atoms:
        for j in metal_ligand_dict_in_cutoff_not_CONECT[i]:
            if not any(topo.atom(j).residue == topo.atom(k).residue for k in metal_ligand_dict_by_CONECT[i]):
                metal_ligand_dict_CORRECTED_in_cutoff_not_CONECT[i].append(j)
    
    
    # compare the cutoff and CONECT data and make conclusions
    for i in metal_atoms:
        
        if metal_ligand_dict_CORRECTED_in_cutoff_not_CONECT[i] and not metal_ligand_dict_in_CONECT_not_cutoff[i]:
            
            #extra_in_cutoff_count_by_atom += 1
            
            #if not extra_in_cutoff_for_file_analyzed:
            #    extra_in_cutoff_count_by_file += 1
            #    extra_in_cutoff_file_list.append(file)
            #    extra_in_cutoff_for_file_analyzed = True
        
            #for j in metal_ligand_dict_CORRECTED_in_cutoff_not_CONECT[i]:
            #    extra_in_cutoff_residue_list.append(topo.atom(j))   
            
            numbers_of_extra_cutoff_ligands.append(len(metal_ligand_dict_CORRECTED_in_cutoff_not_CONECT[i])) 
    
    # update the results dictionary (ok_file_count and error_file_count have already been done)

    #dictionary_of_process_counts['files_with_metal_count'] = files_with_metal_count
    #dictionary_of_process_counts['metal_count'] = metal_count
    #dictionary_of_process_counts['extra_in_cutoff_count_by_atom'] = extra_in_cutoff_count_by_atom
    #dictionary_of_process_counts['extra_in_cutoff_count_by_file'] = extra_in_cutoff_count_by_file
    #dictionary_of_process_counts['extra_in_cutoff_file_list'] = extra_in_cutoff_file_list
    #dictionary_of_process_counts['extra_in_cutoff_residue_list'] = extra_in_cutoff_residue_list
                                                          
    print(file)
    return numbers_of_extra_cutoff_ligands                                                      
                     
#def ligand_scanner_all_database(pdbpath):
    
    #dictionary_of_database_results = {}
    #dictionary_of_database_results['ok_file_count'] = 0
    #dictionary_of_database_results['error_file_count'] = 0
    #dictionary_of_database_results['files_with_metal_count'] = 0
    #dictionary_of_database_results['metal_count'] = 0
    #dictionary_of_database_results['extra_in_cutoff_count_by_atom'] = 0
    #ictionary_of_database_results['extra_in_cutoff_count_by_file'] = 0
    #dictionary_of_database_results['extra_in_cutoff_file_list'] = []
    #dictionary_of_database_results['extra_in_cutoff_residue_list'] = []
    

#    for dictionary_of_process_counts in pool.map(ligand_scanner, glob.iglob(pdbpath)):
        
#        for i in dictionary_of_process_counts:
#            dictionary_of_database_results[i] += dictionary_of_process_counts[i]
    
#    return dictionary_of_database_results
    
# Multiprocess set-up
if __name__ == '__main__':
    
    pool = mp.Pool(processes = ppn)
    results = pool.map(ligand_scanner, glob.iglob(pdbpath))
    

# Write results
f = open('ligand_scanner_results.txt', 'w')
f.write(str(results))
f.close()
