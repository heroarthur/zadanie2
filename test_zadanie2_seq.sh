#!/bin/bash
# FILES="/genomy/"

declare -i index=0

for f in genomy/genom_*
do
  mpirun -n 13 ./zadanie2 $index 0 genomy/genom queries out > f1
  ./seq genomy/genom_$index > f2
  echo "file genom_$index"
  diff f1 f2
  rm f1
  rm f2
  index=$(( index + 1 ))
done