

f1 = open("f1", "r")


f2 = open("f2", "r")

for (x1, x2) in zip(f1.readline().split(), f2.readline().split()):
    if (x1 != x2):
        print("roznia sie ", x1, x2)