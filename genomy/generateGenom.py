import random
import sys
# tr -dc ACGT </dev/urandom | head -c 38472947 > genom_2

fileName = sys.argv[1]
genomeSize = sys.argv[2]


f = open(fileName, 'w+')

genome_text = ''.join([random.choice(['A', 'G', 'T', 'C']) for i in range(int(genomeSize))])

f.write(genome_text)
f.close()




