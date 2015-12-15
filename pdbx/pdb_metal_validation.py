from __future__ import print_function
import gzip
import glob
import multiprocessing as mp
from simtk.openmm import *
from simtk.openmm.app import *
from simtk.unit import *
from ctypes import c_int

ppn = 32
metals = ['Al','As','Ba','Ca','Cd','Ce','Co','Cs','Cu','Dy','Fe','Gd','Hg','Ho','In','Ir','K','Li','Mg',
'Mn','Mo','Na','Ni','Pb','Pd','Pt','Rb','Rh','Ru','Sm','Sr','Te','Tl','V','W','Yb','Zn']

pdblist_file = open('pdblist_conect.txt')
pdblist = [line[:-1] for line in pdblist_file]

def file_reader(file):
    global progress_counter
    with lock:
        progress_counter.value += 1
    print(progress_counter.value)
    
    file_open = gzip.open(file)
    file_reader_results = []
    try:
        pdb = PDBFile(file_open)
    except:
        return file_reader_results
    top = pdb.topology
       
    bonds = [bond for bond in top.bonds()]
    
    for bond in bonds:
        if bond[0].element != None and bond[1].element != None:
            if bond[0].element.symbol in metals or bond[1].element.symbol in metals:
                file_reader_results.append([bond[0].residue.name, bond[0].name, bond[1].residue.name, bond[1].name])   
    
          
    return file_reader_results    
    
    
def database_analyzer(pdblist):
    
    database_analyzer_results = []
    
    for file_reader_results in pool.map(file_reader, pdblist):
        database_analyzer_results.append(file_reader_results)
        
    return database_analyzer_results
        
if __name__ == "__main__":
    lock = mp.Lock()
    progress_counter = mp.Value(c_int)
    pool = mp.Pool(processes = ppn)
    database_analyzer_results = database_analyzer(pdblist) 
    
    # write results
    with open('pdb_metal_validation_results.txt', 'w') as f:
        f.write(str(database_analyzer_results))
                
    
