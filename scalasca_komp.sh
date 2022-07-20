#!/bin/bash
/opt/scalasca/bin/scalasca -instrument --memory --openmp --thread=omp --preprocess mpic++ -fopenmp -O0 --std=c++2a zadanie2.cpp data_source.h data_source.cpp -o zadanie2

