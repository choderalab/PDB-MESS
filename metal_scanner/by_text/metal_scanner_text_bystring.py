from __future__ import print_function
import multiprocessing as mp
import glob
import gzip

pdbpath = '/cbio/jclab/share/pdb/*/*.ent.gz'
ppn = 32
metal_name = 'ZN'

def metal_scanner(file):
    
    contains_metal = False
    pdbfile = gzip.open(file)
    
    for line in pdbfile:
        fields = line.split()
        
        if metal_name in fields:
            contains_metal = True
    
    print(contains_metal)        
    return contains_metal
    
def database_analyzer(pdbpath):
    
    files_w_metal = 0
    files_wo_metal = 0
    
    for contains_metal in pool.map(metal_scanner, glob.iglob(pdbpath)):
        if contains_metal == True:
            files_w_metal += 1
        elif contains_metal == False:
            files_wo_metal += 1    
            
    return files_w_metal, files_wo_metal        
            
# Multiprocess set-up
if __name__ == '__main__':
    
    pool = mp.Pool(processes = ppn)
    files_w_metal, files_wo_metal = database_analyzer(pdbpath)       
    
                    
# write results
with open('metal_scanner_results.txt', 'w') as f:
    f.write('All files analyzed:\n')  
    f.write(str(files_w_metal + files_wo_metal))
    f.write('\n')
    f.write('All files in the database:\n')
    f.write(str(len(glob.glob(pdbpath))))
    f.write('\n')
    f.write('Files containing %s: \n' % metal_name)
    f.write(str(files_w_metal))
    f.write('\n')
                   
