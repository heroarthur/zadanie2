import sys
import re
# python testQueries.py 3 2 genomy/genom genomy/queries out

n = int(sys.argv[1])
m = int(sys.argv[2])
genome_in = sys.argv[3]
queries_in = sys.argv[4]
queries_out = sys.argv[5]

queries = []
results = [[] for i in range(n)]

queriesFile = open(queries_in, "r")
for q in queriesFile:
    queries.append(q.replace('\n', ''))
queriesFile.close()

for g in range(n):
    genomFile = open(genome_in + "_" + str(g), "r")
    genome = genomFile.readline()
    for q in range(m):
        results[g].append(len([index for index in range(len(genome)) if genome.startswith(queries[q], index)]))

f = open(queries_out, "w")

for q in range(m):
    s = ""
    for g in range(n):
        s = s + str(results[g][q])
        s += " "
    s = s[:-1]
    s += '\n'
    f.write(s)

f.close()
