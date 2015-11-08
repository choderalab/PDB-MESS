from __future__ import division
import multiprocessing as mp
import gzip
from collections import Counter

t = open('znconectmissing')
files = t.readlines()
files = [x[:-1] for x in files]
processes = 4
allprops = []
overone = []

def znconect(i):
    
#    global allprops

    for file in i:
        pdbfile = gzip.open(file)
        listofzn = []
        listofznconect = []

        for line in pdbfile:
            fields = line.split()

            if fields[0] == "HETATM" and fields[2] == "ZN" and fields[1] not in listofzn:
                listofzn.append(fields[1])
            elif fields[0] == "CONECT":
                for zn in listofzn:
                    if fields[1] == zn:
                        listofznconect.append(zn)
    
        print(file)
        print(len(listofznconect))
        print(len(listofzn))
    
        allprops.append(len(listofznconect)/len(listofzn))
        if len(listofznconect)/len(listofzn) >= 1:
            overone.append(file)

# multiprocessing
#if __name__ == '__main__':
#    manager = mp.Manager()
#    allprops = manager.list()
#    pool = mp.Pool(processes)
#    pool.map(znconect, files)

znconect(files)

#with open('missinganalresults', 'w') as f:
#    f.write(str(Counter(allprops)))

with open('missinganalresults_overone', 'w') as g:
    for k in overone:
        g.write(str(k))
        g.write("\n")

