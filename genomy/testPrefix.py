#!/bin/bash
mpirun -n 4 ./testMPIprefix 1 0 genom queries out 0 5


# CAACCTTGCGACAGGGCGGG
#         CGACAGGGCTokeanos-login2 zadanie2/genomy> cat testPrefix.py 
import random
import sys
import glob
import subprocess



directory = sys.argv[1]
testCount = int(sys.argv[2])


genomFiles = glob.glob(directory + "/genom_*")

for i in range(testCount):
    genomFile_index = random.randint(0, 17)
    genomFile = "genom_" + str(genomFile_index)
    print(genomFile)
    f = open(genomFile, "r")
    genom = f.readline()
    start = random.randint(0, len(genom)-2)
    end = random.randint(start+1, len(genom)-1)
    f.close()

    bashCommand = f'mpirun -n 4 ./testMPIprefix {genomFile_index} 0 genom queries out {start} {end}'
    process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
    # output, error = process.communicate()
    # prefixMPI = str(output)[2:-1]
    # prefixPython = genom[start:end+1]
    # if (prefixMPI != prefixPython):
    #     print(genomFile)
    #     print("zakresy ", start, end)
    #     break

