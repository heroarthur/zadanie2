#!/bin/bash
rm -dr scorep_zadanie2_*
/opt/scalasca/bin/scalasca -analyze mpirun -n 8 ./zadanie2 17 0 genomy/genom queries out > f
