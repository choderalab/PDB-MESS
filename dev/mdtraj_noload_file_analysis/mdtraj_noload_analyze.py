from __future__ import print_function
import multiprocessing as mp
import glob

pdbpath = 'files/*.ent'
ppn = 6

def file_scanner(file):
    pdbfile = open(file)
    
    lines = [line for line in pdbfile]
    
    conect_lines = []
    conect_atoms = []
    conect_atoms_dictionary = {}
    
    file_scanner_results = [0, 0, 0]
    
    for i in range(len(lines)):
        fields = lines[i].split()
        
        if fields[0] == 'CONECT':
            conect_lines.append(i)
            conect_atoms.append(fields[1])
            
    conect_atoms = list(set(conect_atoms)) 
    
    for conect_atom in conect_atoms:
        conect_atoms_dictionary[conect_atom] = []      
        
    file_scanner_results[2] = len(conect_atoms)     
            
    for line in lines:
        fields = line.split()
        
        if fields[0] == 'ATOM' or fields[0] == 'HETATM':
            for conect_atom in conect_atoms:
                if fields[1] == str(conect_atom):
                    conect_atoms_dictionary[conect_atom].append((fields[2], fields[3]))        
    
    for conect_atom in conect_atoms_dictionary:
        conect_atoms_dictionary[conect_atom] = list(set(conect_atoms_dictionary[conect_atom]))
        
        if len(conect_atoms_dictionary[conect_atom]) != 1:
            file_scanner_results[0] = 1
            file_scanner_results[1] += 1
            
    return file_scanner_results

def database_analyzer(pdbpath):
    
    problem_file_counter = 0
    problem_atom_tuples = []
    
    for file_scanner_results in pool.map(file_scanner, glob.iglob(pdbpath)):
        if file_scanner_results[0] == 1:
            problem_file_counter += 1
            problem_atom_tuples.append((file_scanner_results[1], file_scanner_results[2]))
            
    database_analyzer_results = [problem_file_counter, problem_atom_tuples]        
            
    return database_analyzer_results 
            
        
if __name__ == '__main__':
    pool = mp.Pool(processes = ppn)        
        
    database_analyzer_results = database_analyzer(pdbpath)
    problem_file_counter = database_analyzer_results[0]
    problem_atom_tuples = database_analyzer_results[1]
    
with open('mdtraj_noload_analyze_results.txt') as f:
    f.write('problem_file_counter\n')
    f.write(str(problem_file_counter))
    f.write('\n')
    f.write('problem_atom_tuples\n')
    f.write(problem_atom_tuples)        
