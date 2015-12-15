from __future__ import print_function
import gzip
import glob
import multiprocessing as mp

database_path = '/cbio/jclab/share/pdbx/*/*.cif.gz'
ppn = 32


def file_reader(file):
    file_open = gzip.open(file)
    file_reader_results = []

    
    for line in file_open:
        fields = line.split()
        
        if len(fields) == 33 and (fields[1] == 'covale' or fields[1] == 'modres'):
            file_reader_results.append([fields[4], fields[6], fields[12], fields[14]])
    
    
    
    print(file_reader_results)        
    return file_reader_results    
    
    
def database_analyzer(database_path):
    
    database_analyzer_results = []
    
    for file_reader_results in pool.map(file_reader, glob.iglob(database_path)):
        database_analyzer_results.append(file_reader_results)
        
    return database_analyzer_results
        
if __name__ == "__main__":
    pool = mp.Pool(processes = ppn)
    database_analyzer_results = database_analyzer(database_path) 
    
    # write results
    with open('pdbx_conect_reader_covale_results.txt', 'w') as f:
        f.write(str(database_analyzer_results))
                
    
