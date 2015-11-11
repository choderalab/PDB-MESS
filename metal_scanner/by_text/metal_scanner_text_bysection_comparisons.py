from __future__ import print_function
import multiprocessing as mp
import glob
import gzip

pdbpath = '/cbio/jclab/share/pdb/*/*.ent.gz'
ppn = 32
metal_name = 'ZN'

def metal_scanner(file):
    
    metal_scanner_data = [0, 0, 0, file]
    pdbfile = gzip.open(file)
    
    for line in pdbfile:
        fields = line.split()
            
        if fields[0] == 'HET' and fields[1] == metal_name:
            metal_scanner_data[0] = 1
        if fields[0] == 'HETATM' and fields[2] == metal_name:
            metal_scanner_data[1] = 1    
        if 'ZN' in fields and 'REMARK' not in fields:
            metal_scanner_data[2] = 1    
    
    print(metal_scanner_data)        
    return metal_scanner_data
    
def database_analyzer(pdbpath):
    
# Questions to answer: Do the HET and HETATM sets all contain within the 'ZN' in line sets, is HETATM all contained within HET, give me all the 'ZN' in line that do not have HET or HETAMT
# Also give me all HETs that do not have HETATMs 
# So: are there any [x,x,0] where x is 1    
    
    # HET and HETATM not within 'ZN'-line , HETATM not in HET, HET without HETATM - file path, ZN-line without HET and HETATM - file path, 
    database_analyzer_output = [0, 0, [], []]
    
    for metal_scanner_data in pool.map(metal_scanner, glob.iglob(pdbpath)):
        
        if metal_scanner_data[0] == 1 or metal_scanner_data[1] == 1:
            if metal_scanner_data[2] == 0:
                database_analyzer_output[0] += 1
                
        if metal_scanner_data[1] == 1:
            if metal_scanner_data[0] == 0:
                database_analyzer_output[1] += 1
                
        if metal_scanner_data[0] == 1 and metal_scanner_data[1] == 0:
            database_analyzer_output[2].append(metal_scanner_data[3])
            
        if metal_scanner_data[2] == 1 and metal_scanner_data[0] == 0 and metal_scanner_data[1] == 0:
            database_analyzer_output[3].append(metal_scanner_data[3])          
                        
    return database_analyzer_output     
            
# Multiprocess set-up
if __name__ == '__main__':
    
    pool = mp.Pool(processes = ppn)
    database_analyzer_output = database_analyzer(pdbpath)       
    
                    
# write results
with open('metal_scanner_results_by_text_bysection_comparison.txt', 'w') as f:
    f.write("metal_name\n")
    f.write(metal_name)
    f.write("\n")
    f.write('All files in the database:\n')
    f.write(str(len(glob.glob(pdbpath))))
    f.write('\n')
    f.write('HET and HETATM not contained within ZN-line:\n')
    f.write(str(database_analyzer_output[0]))
    f.write('\n')
    f.write('HETATM not contained within HET:\n')
    f.write(str(database_analyzer_output[1]))
    f.write('\n')
    f.write('HET without HETATM files:\n')
    f.write('Number:\n')
    f.write(str(len(database_analyzer_output[2])))
    f.write('\n')
    for i in database_analyzer_output[2]:
        f.write(str(i))
        f.write('\n')
    f.write('ZN-line without HET or HETATM:\n')
    f.write('Number:\n')
    f.write(str(len(database_analyzer_output[3])))
    f.write('\n')
    for i in database_analyzer_output[3]:    
        f.write(str(i))
        f.write('\n')
                   
