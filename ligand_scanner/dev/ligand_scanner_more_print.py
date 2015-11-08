from __future__ import print_function
import mdtraj as md
import glob
import multiprocessing as mp
from ctypes import c_int

pdbpath = '/cbio/jclab/share/pdb/*/*.ent.gz'
cutoff = 0.3
metal_name = 'ZN'
ppn = 32

def ligand_scanner(file):
    
    global ok_file_count, error_file_count, metal_count, files_with_metal_count
    global more_than_CONECT_count_by_atom, less_than_CONECT_count_by_atom, equal_CONECT_count_by_atom
    global more_than_CONECT_count_by_file, less_than_CONECT_count_by_file, equal_CONECT_count_by_file
    
    metal_ligand_dict_by_cutoff = {}
    metal_ligand_dict_by_CONECT = {}
    equal_for_file_analyzed = False
    more_for_file_analyzed = False
    less_for_file_analyzed = False
    
    # print file path into terminal
    #print(file)
    
    try:
        traj = md.load_pdb(file)
        with lock:
            ok_file_count.value += 1
    except:
        with lock:
            error_file_count.value += 1
        return None

    
    # Select atoms and pairs
    topo = traj.topology
    metal_select_name = 'name %s' % metal_name
    metal_atoms = topo.select(metal_select_name)
    metal_all_pairs = topo.select_pairs(metal_select_name, 'all')
    
    # No metal - skip, metal - add to count 
    if not metal_atoms.size:
        return None
    
    files_with_metal_count.value += 1
    
    # Prep dictionary with metal indices as keys
    for i in metal_atoms:
        metal_ligand_dict_by_cutoff[i] = []
        metal_ligand_dict_by_CONECT[i] = []
        with lock:
            metal_count.value += 1
        
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
                
    # compare the cutoff to CONECT 
    for i in metal_ligand_dict_by_cutoff:
        
        #if len(metal_ligand_dict_by_cutoff[i]) == len(metal_ligand_dict_by_CONECT[i]):
        #    with lock:
    #            equal_CONECT_count_by_atom.value += 1
#            if not equal_for_file_analyzed:
#               with lock:
            #       equal_CONECT_count_by_file.value += 1  
             #  equal_for_file_analyzed = True 
        
        #if len(metal_ligand_dict_by_cutoff[i]) > len(metal_ligand_dict_by_CONECT[i]):
            with lock:
                more_than_CONECT_count_by_atom.value += 1  
            if not more_for_file_analyzed:
               with lock:
                   more_than_CONECT_count_by_file.value += 1 
               print(file)
               more_for_file_analyzed = True
        
        #if len(metal_ligand_dict_by_cutoff[i]) < len(metal_ligand_dict_by_CONECT[i]):
            with lock:
                less_than_CONECT_count_by_atom.value += 1     
            if not less_for_file_analyzed:
               with lock:
                   less_than_CONECT_count_by_file.value += 1 
               less_for_file_analyzed = True                 
               
               
# Multiprocess
if __name__ == '__main__':
    manager = mp.Manager()
    lock = mp.Lock()
    ok_file_count = mp.Value(c_int)
    error_file_count = mp.Value(c_int)
    files_with_metal_count = mp.Value(c_int)
    metal_count = mp.Value(c_int)
    more_than_CONECT_count_by_atom = mp.Value(c_int)
    less_than_CONECT_count_by_atom = mp.Value(c_int)
    equal_CONECT_count_by_atom = mp.Value(c_int)
    more_than_CONECT_count_by_file = mp.Value(c_int)
    less_than_CONECT_count_by_file = mp.Value(c_int)
    equal_CONECT_count_by_file = mp.Value(c_int)

    pool = mp.Pool(processes = ppn)
    pool.map(ligand_scanner, glob.iglob(pdbpath))


# Write results
#f = open('ligand_scanner_results.txt', 'w')
#f.write("ok_file_count: " + str(ok_file_count.value) + "\n")
#f.write("error_file_count: " + str(error_file_count.value) + "\n")
#f.write("files_with_metal_count: " + str(files_with_metal_count.value) + "\n")
#f.write("metal_count: " + str(metal_count.value) + "\n")
#f.write("more_than_CONECT_count_by_atom: " + str(more_than_CONECT_count_by_atom.value) + "\n")
#f.write("less_than_CONECT_count_by_atom: " + str(less_than_CONECT_count_by_atom.value) + "\n")
#f.write("equal_CONECT_count_by_atom: " + str(equal_CONECT_count_by_atom.value) + "\n")
#f.write("more_than_CONECT_count_by_file: " + str(more_than_CONECT_count_by_file.value) + "\n")
#f.write("less_than_CONECT_count_by_file: " + str(less_than_CONECT_count_by_file.value) + "\n")
#f.write("equal_CONECT_count_by_file: " + str(equal_CONECT_count_by_file.value) + "\n")
#f.close()
