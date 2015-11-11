from __future__ import print_function
import multiprocessing as mp
import glob
import gzip

pdbpath = '/cbio/jclab/share/pdb/*/*.ent.gz'
ppn = 32
metal_name = 'ZN'

def metal_scanner(file):
    
    metal_scanner_data = [0,0]
    pdbfile = gzip.open(file)
    
    for line in pdbfile:
        fields = line.split()
        
        if fields[0] == 'HET' and fields[1] == metal_name:
            metal_scanner_data[0] = 1
        if fields[0] == 'HETATM' and fields[2] == metal_name:
            metal_scanner_data[1] = 1    
    
    print(metal_scanner_data)        
    return metal_scanner_data
    
def database_analyzer(pdbpath):
    
    files_w_HET = 0
    files_w_HETATM = 0
    
    for metal_scanner_data in pool.map(metal_scanner, glob.iglob(pdbpath)):
        if metal_scanner_data[0] == 1:
            files_w_HET += 1
        if metal_scanner_data[1] == 1:
            files_w_HETATM += 1       
            
    return files_w_HET, files_w_HETATM       
            
# Multiprocess set-up
if __name__ == '__main__':
    
    pool = mp.Pool(processes = ppn)
    files_w_HET, files_w_HETATM = database_analyzer(pdbpath)       
    
                    
# write results
with open('metal_scanner_results_by_text_bysection.txt', 'w') as f:
    f.write("metal_name\n")
    f.write(metal_name)
    f.write("\n")
    f.write('All files in the database:\n')
    f.write(str(len(glob.glob(pdbpath))))
    f.write('\n')
    f.write("Files with metal in HET:\n")
    f.write(str(files_w_HET))
    f.write('\n')
    f.write("Files with metal in HETATM:\n")
    f.write(str(files_w_HETATM))
    f.close()
                   
