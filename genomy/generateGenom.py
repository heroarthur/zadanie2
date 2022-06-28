import random
import sys


fileName = sys.argv[1]
genomeSize = sys.argv[2]


f = open(fileName, 'w+')

for i in range(int(genomeSize)):
    f.write(random.choice(['A', 'G', 'T', 'C']))

f.close()