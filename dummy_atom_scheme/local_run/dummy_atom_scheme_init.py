from __future__ import print_function
import mdtraj as md
import glob
import multiprocessing as mp

pdblist_file = open('pdblist.txt')
pdblist = [line[:-1] for line in pdblist_file]
cutoff = 0.3
metal_name = 'ZN'
ppn = 8


def ligand_scanner(file):
    
    print(file)
    
    ok_file_count, error_file_count, metal_count, files_with_metal_count = 0, 0, 0, 0
    metal_ligand_dict_by_cutoff = {}
    metal_ligand_dict_by_CONECT = {}
    ligand_numbers = []
    CONECT_atoms = []
    CUTOFF_atoms = []
    double_atoms = []
    CONECT_residues = []
    CUTOFF_residues = []
    metal_ligand_dict_cutoff_analysis = {}
    metal_ligand_cutoff_analysis = []
    
    dictionary_of_process_counts = {}
    
    # initialize the results dictionary (dictionary_of_process_counts)        
    dictionary_of_process_counts['ok_file_count'] = ok_file_count
    dictionary_of_process_counts['error_file_count'] = error_file_count
    dictionary_of_process_counts['files_with_metal_count'] = files_with_metal_count
    dictionary_of_process_counts['metal_count'] = metal_count
    dictionary_of_process_counts['ligand_numbers'] = ligand_numbers
    dictionary_of_process_counts['CONECT_atoms'] = CONECT_atoms
    dictionary_of_process_counts['CUTOFF_atoms'] = CUTOFF_atoms
    dictionary_of_process_counts['double_atoms'] = double_atoms
    dictionary_of_process_counts['CONECT_residues'] = CONECT_residues
    dictionary_of_process_counts['CUTOFF_residues'] = CUTOFF_residues
    dictionary_of_process_counts['metal_ligand_cutoff_analysis'] = metal_ligand_cutoff_analysis
    
    
    
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
    metal_all_pairs = topo.select_pairs(metal_select_name, 'symbol O or symbol N or symbol S or symbol Cl or symbol F')
    
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
                
                
                
    # Further analysis on CUTOFF - 1) Find: HOH-O, CYS-S, HIS-N, ASP-O, GLU-O (but only one for those) - so count as 'residues' rather than atoms (to correct for the ASP and GLU)
    # 2) See if there are any extras in the cutoff of elements O, Cl or F
    
    for i in metal_atoms:
        
        metal_ligand_dict_cutoff_analysis[i] = [[], 0, 0]
        
        for j in metal_ligand_dict_by_cutoff[i]:
            metal_ligand_dict_cutoff_analysis[i][0].append((topo.atom(j).residue, topo.atom(j).residue.chain))
        metal_ligand_dict_cutoff_analysis[i][0] = list(set(metal_ligand_dict_cutoff_analysis[i][0]))    
        
        for k in metal_ligand_dict_cutoff_analysis[i][0]:
            if k[0].name == 'ASP' or k[0].name == 'GLU' or k[0].name == 'HOH' or k[0].name == 'HIS' or k[0].name == 'CYS':
                metal_ligand_dict_cutoff_analysis[i][1] += 1
            else:
                metal_ligand_dict_cutoff_analysis[i][2] += 1    
                
    for i in metal_atoms:
        metal_ligand_cutoff_analysis.append(metal_ligand_dict_cutoff_analysis[i])
        
                    
    # Compare the performance of CONECT and CUTOFF
    for i in metal_atoms:
        ligand_numbers.append((file, i, len(metal_ligand_dict_by_CONECT[i]), len(metal_ligand_dict_by_cutoff[i])))
        
        for j in metal_ligand_dict_by_CONECT[i]:
            CONECT_atoms.append((file, i, str(topo.atom(j)), topo.atom(j).element))
            for k in (k for k in metal_ligand_dict_by_CONECT[i] if k != j):
                if topo.atom(j).residue == topo.atom(k).residue:
                    double_atoms.append((str(topo.atom(j)), str(topo.atom(k))))
            CONECT_residues.append((file, i, str(topo.atom(j).residue.name)))        
            
        
        for j in metal_ligand_dict_by_cutoff[i]:
            CUTOFF_atoms.append((file, i, str(topo.atom(j)), topo.atom(j).element))
            CUTOFF_residues.append((file, i, str(topo.atom(j).residue.name)))
        

    # update the results dictionary (ok_file_count and error_file_count have already been done)

    dictionary_of_process_counts['ok_file_count'] = ok_file_count
    dictionary_of_process_counts['error_file_count'] = error_file_count
    dictionary_of_process_counts['files_with_metal_count'] = files_with_metal_count
    dictionary_of_process_counts['metal_count'] = metal_count
    dictionary_of_process_counts['ligand_numbers'] = ligand_numbers
    dictionary_of_process_counts['CONECT_atoms'] = CONECT_atoms
    dictionary_of_process_counts['CUTOFF_atoms'] = CUTOFF_atoms
    dictionary_of_process_counts['double_atoms'] = double_atoms
    dictionary_of_process_counts['CONECT_residues'] = CONECT_residues
    dictionary_of_process_counts['CUTOFF_residues'] = CUTOFF_residues
    dictionary_of_process_counts['metal_ligand_cutoff_analysis'] = metal_ligand_cutoff_analysis
                
                                                          
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
    dictionary_of_database_results['double_atoms'] = []
    dictionary_of_database_results['CONECT_residues'] = []
    dictionary_of_database_results['CUTOFF_residues'] = []
    dictionary_of_database_results['metal_ligand_cutoff_analysis'] = []
     
    
    for dictionary_of_process_counts in pool.map(ligand_scanner, pdblist):
        for i in dictionary_of_process_counts:
            dictionary_of_database_results[i] += dictionary_of_process_counts[i]
            
    return dictionary_of_database_results
    
# Multiprocess set-up
if __name__ == '__main__':
    
    pool = mp.Pool(processes = ppn)
    dictionary_of_database_results = ligand_scanner_all_database(pdblist)
    

# Write results
f = open('ligand_scanner_towards_dummy_scheme.txt', 'w')
f.write(str(dictionary_of_database_results))
f.close()
