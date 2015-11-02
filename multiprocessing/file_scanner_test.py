import mdtraj as md
import glob
import multiprocessing as mp
from ctypes import c_int

pdbpath = '/cbio/jclab/share/pdb/*/*.ent.gz'

def file_reader(file):
    
    global ok_file_count
    global error_file_count
    
    try:
        traj = md.load_pdb(file)
        with lock:
            ok_file_count.value += 1
    except:
        with lock:
            error_file_count.value += 1
    
    print(file)
    print(ok_file_count.value)        
        
    
if __name__ == '__main__':
    manager = mp.Manager()
    lock = mp.Lock()
    ok_file_count = mp.Value(c_int)
    error_file_count = mp.Value(c_int)
    pool = mp.Pool(processes = 25)
    pool.map(file_reader, glob.iglob(pdbpath))        
        
f = open('file_reader_test_results.txt', 'w')
f.write("ok_file_count: " + str(ok_file_count.value) + "\n")
f.write("error_file_count: " + str(error_count_file.value))
f.close()        
        
