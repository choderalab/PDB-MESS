from __future__ import print_function
import gzip
import glob
import multiprocessing as mp
from simtk.openmm import *
from simtk.openmm.app import *
from simtk.unit import *
from ctypes import c_int

database_path = '/cbio/jclab/share/pdb/*/*.ent.gz'
ppn = 32


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
    
    for atom in top.atoms():
        if atom.element == None:
            if not file_reader_results:
                file_reader_results.append(file)
            file_reader_results.append(atom.residue.name)
    
    return file_reader_results    
    
    
def database_analyzer(database_path):
    
    database_analyzer_results = []
    
    for file_reader_results in pool.map(file_reader, glob.iglob(database_path)):
        database_analyzer_results.append(file_reader_results)
        
    return database_analyzer_results
        
if __name__ == "__main__":
    lock = mp.Lock()
    progress_counter = mp.Value(c_int)
    pool = mp.Pool(processes = ppn)
    database_analyzer_results = database_analyzer(database_path) 
    
    # write results
    with open('PDB_reader_elementcheck_results.txt', 'w') as f:
        f.write(str(database_analyzer_results))
                
    
