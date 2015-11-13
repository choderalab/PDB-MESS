from __future__ import print_function
import multiprocessing as mp
import glob
import mdtraj as md

pdbpath = '/cbio/jclab/share/pdb/*/*.ent.gz'
ppn = 32
metal_name = 'ZN'

def metal_scanner(file):
    
    contains_metal = [False, file]
    try:
        traj = md.load_pdb(file)
    except:
        return [None, file]    
    
    if [atom for atom in traj.top.atoms if atom.name == 'ZN']:
      contains_metal[0] = True
    
    print(contains_metal)        
    return contains_metal
    
def database_analyzer(pdbpath):
    
    database_analyzer_results = [0, []]
    #files_wo_metal = 0
    
    for contains_metal in pool.map(metal_scanner, glob.iglob(pdbpath)):
        if contains_metal[0] == True:
            database_analyzer_results[0] += 1
            database_analyzer_results[1].append(contains_metal[1])
    #    elif contains_metal == False:
        #    files_wo_metal += 1    
            
    return database_analyzer_results       
            
# Multiprocess set-up
if __name__ == '__main__':
    
    pool = mp.Pool(processes = ppn)
    database_analyzer_results = database_analyzer(pdbpath)       
    
                    
# write results
with open('metal_scanner_by_mdtraj_results_FILES.txt', 'w') as f:
    #f.write('All files analyzed:\n')  
    #f.write(str(files_w_metal + files_wo_metal))
    #f.write('\n')
    #f.write('All files in the database:\n')
    #f.write(str(len(glob.glob(pdbpath))))
    #f.write('\n')
    f.write('Files containing %s: \n' % metal_name)
    f.write(str(database_analyzer_results[0]))
    f.write('\n')
    f.write('File paths:\n')
    for line in database_analyzer_results[1]:
        f.write(str(line))
        f.write('\n')
        
                   
