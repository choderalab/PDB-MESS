from __future__ import print_function
import multiprocessing as mp
import glob
import gzip

metal_name = 'ZN'
ppn = 32
pdbpath = '/cbio/jclab/share/pdb/*/*.ent.gz'

def file_scanner(file):
    
    pdbfile = gzip.open(file)
    metal_counter = 0
    
    
    for line in pdbfile:
        fields = line.split()
        
        if any(field == metal_name for field in fields):
            metal_counter = 1
            
    print(metal_counter)        
    return metal_counter        

def database_scanner(pdbpath):
    
    metal_counter_global = 0
        
    for metal_counter in pool.map(file_scanner, glob.iglob(pdbpath)):
        if metal_counter == 1:
            metal_counter_global += 1
            
    return metal_counter_global
    
if __name__ == '__main__':
    
    pool = mp.Pool(processes = ppn)
    metal_counter_global = database_scanner(pdbpath)
    
f = open('database_scanner_results.txt', 'w')
f.write(str(metal_counter_global))
f.close()
    
                
                        
