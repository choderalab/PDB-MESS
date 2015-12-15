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
    dictionary_of_process_counts = {}
    ligand_numbers = []
    five_coord_distances = []
    five_coord_conect = []
    
    # initialize the results dictionary (dictionary_of_process_counts)        
    dictionary_of_process_counts['ok_file_count'] = ok_file_count
    dictionary_of_process_counts['error_file_count'] = error_file_count
    #dictionary_of_process_counts['files_with_metal_count'] = files_with_metal_count
    #dictionary_of_process_counts['metal_count'] = metal_count
    dictionary_of_process_counts['ligand_numbers'] = ligand_numbers
    dictionary_of_process_counts['five_coord_distances'] = five_coord_distances
    dictionary_of_process_counts['five_coord_conect'] = five_coord_conect
    
    
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
    hoh_ligand_atoms = [atom.index for atom in topo.atoms if atom.name == 'O' and atom.residue.name == 'HOH']
    other_ligand_atoms = [atom.index for atom in topo.atoms if (atom.name =="CL" and atom.residue.name == "CL") or (atom.name == 'F' and atom.residue.name == 'F') or (tuple(atom.element)[2] == "O" and atom.residue.is_protein == False)]
    
    # no water and others for now
    ligand_atoms = cys_ligand_atoms + his_ligand_atoms + glu_ligand_atoms + asp_ligand_atoms
    
    metal_all_pairs = topo.select_pairs(metal_atoms, ligand_atoms)
    
    # No metal - skip, metal - add to count 
    #if not metal_atoms.size:
    #    return dictionary_of_process_counts
    
    #files_with_metal_count += 1
    #metal_count += len(metal_atoms)
    
    # Prep dictionary with metal indices as keys
    for i in metal_atoms:
        metal_ligand_dict_by_cutoff[i] = []
        
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
        metal_ligand_dict_by_CONECT[i] = []
    
    
    for i in metal_atoms:
        for bond in topo.bonds:
            if bond[0].index == i:
                metal_ligand_dict_by_CONECT[i].append(bond[1].index)
            elif bond[1].index == i:
                metal_ligand_dict_by_CONECT[i].append(bond[0].index)   
    
    
    # Ligand number analysis
    
    for i in metal_ligand_dict_by_cutoff:
        ligand_analysis = []
        five_coord_distances_per_atom = []
        five_coord_conect_per_atom = [[], []]
        
        for j in metal_ligand_dict_by_cutoff[i]:
            ligand_analysis.append((topo.atom(j).residue, topo.atom(j).residue.chain))
        
        ligand_analysis = list(set(ligand_analysis))       
        
        ligand_numbers.append(len(ligand_analysis))
        
        if len(ligand_analysis) == 5:
            for k in ligand_analysis:
                for j in (j for j in metal_ligand_dict_by_cutoff[i] if k[0] == topo.atom(j).residue):
                    pair = topo.select_pairs([i], [j])
                    distance = md.compute_distances(traj, pair, periodic = False)
                    five_coord_distances_per_atom.append((str(topo.atom(j).residue), topo.atom(j).residue.name, distance))
            
            # compare to conect
           
            five_coord_conect_per_atom[0] = ligand_analysis       
        
            for k in metal_ligand_dict_by_CONECT[i]:
                five_coord_conect_per_atom[1].append(str(topo.atom(k)))      
                    
        five_coord_distances.append(five_coord_distances_per_atom)
        five_coord_conect.append(five_coord_conect_per_atom)            
                

        
    # Update results dictionary
    dictionary_of_process_counts['ok_file_count'] = ok_file_count
    dictionary_of_process_counts['error_file_count'] = error_file_count
    #dictionary_of_process_counts['files_with_metal_count'] = files_with_metal_count
    #dictionary_of_process_counts['metal_count'] = metal_count
    dictionary_of_process_counts['ligand_numbers'] = ligand_numbers
    dictionary_of_process_counts['five_coord_distances'] = five_coord_distances
    dictionary_of_process_counts['five_coord_conect'] = five_coord_conect
    
    return dictionary_of_process_counts
    
    
def ligand_scanner_all_database(pdblist):
    
    dictionary_of_database_results = {}
    dictionary_of_database_results['ok_file_count'] = 0
    dictionary_of_database_results['error_file_count'] = 0
    #dictionary_of_database_results['files_with_metal_count'] = 0
    #dictionary_of_database_results['metal_count'] = 0
    dictionary_of_database_results['ligand_numbers'] = []
    dictionary_of_database_results['five_coord_distances'] = []
    dictionary_of_database_results['five_coord_conect'] = []

    
    for dictionary_of_process_counts in map(ligand_scanner, pdblist):
        for i in dictionary_of_process_counts:
            dictionary_of_database_results[i] += dictionary_of_process_counts[i]
            
    return dictionary_of_database_results
    
# Multiprocess set-up
if __name__ == '__main__':
    
    pool = mp.Pool(processes = ppn)
    dictionary_of_database_results = ligand_scanner_all_database(pdblist)
    

# Write results
f = open('debugging.txt', 'w')
f.write(str(dictionary_of_database_results))
f.close()        
                
                
                
