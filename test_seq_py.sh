#!/bin/bash
# FILES="/genomy/"

for f in genomy/genom_*
do
  python buildIndex.py $f > f1
  ./seq $f > f2
  echo "file $f"
  diff f1 f2
  rm f1
  rm f2
done