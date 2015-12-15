from __future__ import print_function
import gzip
import glob
import multiprocessing as mp

database_path = '/cbio/jclab/share/pdb/*/*.ent.gz'
ppn = 32

def file_reader(file):
    print(file)   
    file_open = gzip.open(file)
    
    for line in file_open:
        
        if line[:6] == 'CONECT':
            print('CONECT')
            return file
    
     
def database_analyzer(database_path):
    
    database_analyzer_results = []
    
    for file_reader_results in pool.map(file_reader, glob.iglob(database_path)):
        database_analyzer_results.append(file_reader_results)
        
    return database_analyzer_results
        
if __name__ == "__main__":
    pool = mp.Pool(processes = ppn)
    database_analyzer_results = database_analyzer(database_path) 
    
    # write results
    with open('pdb_conect_files_results.txt', 'w') as f:
        f.write(str(database_analyzer_results))
                
    
