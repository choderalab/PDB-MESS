from __future__ import print_function
import mdtraj as md
import glob
import multiprocessing as mp
from numpy import mean

pdblist_file = open('pdblist.txt')
pdblist = [line[:-1] for line in pdblist_file]
cutoff = 0.3
metal_name = 'ZN'
ppn = 32


def ligand_scanner(file):
    
    print(file)
    
    ok_file_count, error_file_count, metal_count, files_with_metal_count = 0, 0, 0, 0
    metal_ligand_dict_by_cutoff = {}
    dictionary_of_process_counts = {}
    ligand_numbers = []
    over_4_files = []
    over_four_analysis = []
    
    # initialize the results dictionary (dictionary_of_process_counts)        
    dictionary_of_process_counts['ok_file_count'] = ok_file_count
    dictionary_of_process_counts['error_file_count'] = error_file_count
    #dictionary_of_process_counts['files_with_metal_count'] = files_with_metal_count
    #dictionary_of_process_counts['metal_count'] = metal_count
    dictionary_of_process_counts['ligand_numbers'] = ligand_numbers
    dictionary_of_process_counts['over_four_analysis'] = over_four_analysis
    
    
    
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
    #other_ligand_atoms = [atom.index for atom in topo.atoms if (tuple(atom.element)[2] == 'O' or tuple(atom.element)[2] == 'Cl' or tuple(atom.element)[2] == 'F') and atom.residue.name not in ['CYS', 'HIS', 'GLU', 'ASP', 'HOH']]
    other_ligand_atoms = [atom.index for atom in topo.atoms if (atom.name =="CL" and atom.residue.name == "CL") or (atom.name == 'F' and atom.residue.name == 'F') or (tuple(atom.element)[2] == "O" and atom.residue.is_protein == False)]

    
    ligand_atoms = cys_ligand_atoms + his_ligand_atoms + glu_ligand_atoms + asp_ligand_atoms + hoh_ligand_atoms + other_ligand_atoms
    
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
                
    # Ligand number analysis
    
    for i in metal_ligand_dict_by_cutoff:
        ligand_analysis = []
        ligand_numbers_per_metal = [0, 0, 0]
        over_four_analysis_per_metal = [[], []]
        
        for j in metal_ligand_dict_by_cutoff[i]:
            ligand_analysis.append((topo.atom(j).residue, topo.atom(j).residue.chain))
        
        ligand_analysis = list(set(ligand_analysis))      
        
        for k in ligand_analysis: 
            if k[0].name in ['CYS', 'HIS', 'ASP', 'GLU']:
                ligand_numbers_per_metal[0] += 1
            elif k[0].name == 'HOH':
                ligand_numbers_per_metal[1] += 1
            else:
                ligand_numbers_per_metal[2] += 1        
                
                
        ligand_numbers.append(ligand_numbers_per_metal)
        
        if ligand_numbers_per_metal[0] + ligand_numbers_per_metal[1] + ligand_numbers_per_metal[2] > 4:
            for k in ligand_analysis:
                over_four_analysis_per_metal[0].append(k[0].name)
                for x in [x for x in metal_ligand_dict_by_cutoff[i] if k == (topo.atom(x).residue, topo.atom(x).residue.chain)]:
                    distances = []
                    pair = topo.select_pairs([i], [x])
                    distance = md.compute_distances(traj, pair, periodic = False)
                    distances.append(distance)
                over_four_analysis_per_metal[1].append(float(mean(distance)))
            
        over_four_analysis.append(over_four_analysis_per_metal)    
            
        
        
    # Update results dictionary
    dictionary_of_process_counts['ok_file_count'] = ok_file_count
    dictionary_of_process_counts['error_file_count'] = error_file_count
    #dictionary_of_process_counts['files_with_metal_count'] = files_with_metal_count
    #dictionary_of_process_counts['metal_count'] = metal_count
    dictionary_of_process_counts['ligand_numbers'] = ligand_numbers
    dictionary_of_process_counts['over_four_analysis'] = over_four_analysis
    
    return dictionary_of_process_counts
    
    
def ligand_scanner_all_database(pdblist):
    
    dictionary_of_database_results = {}
    dictionary_of_database_results['ok_file_count'] = 0
    dictionary_of_database_results['error_file_count'] = 0
    #dictionary_of_database_results['files_with_metal_count'] = 0
    #dictionary_of_database_results['metal_count'] = 0
    dictionary_of_database_results['ligand_numbers'] = []
    dictionary_of_database_results['over_four_analysis'] = []
    
    for dictionary_of_process_counts in pool.map(ligand_scanner, pdblist):
        for i in dictionary_of_process_counts:
            dictionary_of_database_results[i] += dictionary_of_process_counts[i]
            
    return dictionary_of_database_results
    
# Multiprocess set-up
if __name__ == '__main__':
    
    pool = mp.Pool(processes = ppn)
    dictionary_of_database_results = ligand_scanner_all_database(pdblist)
    

# Write results
f = open('02Dec_dummy_atom_scheme_v4.txt', 'w')
f.write(str(dictionary_of_database_results))
f.close()        
                
                
                
