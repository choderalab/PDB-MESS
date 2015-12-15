from __future__ import print_function
import gzip
import glob
import multiprocessing as mp
from simtk.openmm import *
from simtk.openmm.app import *
from simtk.unit import *
from ctypes import c_int

ppn = 32

pdblist_file = open('pdbxlist_metalc.txt')
pdblist = [line[:-1] for line in pdblist_file]

def file_reader(file):
    global progress_counter
    with lock:
        progress_counter.value += 1
    print(progress_counter.value)
    
    file_open = gzip.open(file)
    file_reader_results = []
    try:
        pdb = PDBxFile(file_open)
    except:
        return file_reader_results
    top = pdb.topology
    
    for atom in top.atoms():
        if atom.element != None:
            if atom.element.symbol not in ['C', 'N', 'O', 'H', 'S', 'Cl', 'Br', 'F', 'Se']:
                file_reader_results.append(atom.element.symbol)    
    
          
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
    with open('pdbx_metalc_elements.txt', 'w') as f:
        f.write(str(database_analyzer_results))
                
    
