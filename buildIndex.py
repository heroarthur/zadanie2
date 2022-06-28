import numpy as np
from random import randint
from collections import deque
import random
import sys

fileName = sys.argv[1]

f = open(fileName, "r")


genome = deque(list(f.read() + '$'))

l = []

for i in range(len(genome)):
    l.append((''.join(genome), i))
    genome.popleft()


l = list(l)
l.sort()
l = [str(elem[1]) for elem in l]

SA_output = ' '.join(l)
print(SA_output)
