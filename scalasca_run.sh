#!/bin/bash
rm -dr scorep_zadanie2_*
/opt/scalasca/bin/scalasca -analyze mpirun -n 1 ./zadanie2 16 0 genomy/genom queries out > f
