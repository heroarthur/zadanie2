#!/bin/bash

declare -i index=0

for f in genomy/genom_*
do
  mpirun -n 4 ./zadanie2 $index 0 genomy/genom queries out > /dev/null
  echo "test $index"
  index=$(( index + 1 ))
done