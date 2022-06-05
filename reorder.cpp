#include <iostream> 
#include <stdio.h>
#include <algorithm>
#include <stdlib.h>
#include <time.h>
#include <cstdio>
#include <cstring>
#include <vector>

#include <omp.h>
#include <mpi.h>

#ifndef   sortsOpenMP
#define   sortsOpenMP

#include "sortsOpenMP.cpp"

#endif


#include<bits/stdc++.h>

#define root 0


using namespace std;


typedef struct ad {
    int64 rank;
    int64 size;
} ArrData;

MPI_Datatype MPI_ArrData;



void arraysSizes(vector<int64>* A, vector<int64>* A_sizes, vector<int64>* A_offsets,  int rank, int worldSize) {
    ArrData thisArr;
    thisArr.rank = rank;
    thisArr.size = A->size();

    vector<ArrData> buffWithSizes; buffWithSizes.resize(worldSize);

    MPI_Allgather(&thisArr, 1, MPI_ArrData, buffWithSizes.data(), worldSize, MPI_ArrData, MPI_COMM_WORLD);

    sort(buffWithSizes.begin(), buffWithSizes.end(), [](ArrData a, ArrData b) { return a.rank < b.rank; });

    A_offsets->data()[0] = 0;
    A_sizes->data()[0] = buffWithSizes.data()[0].size;
    for (int i = 1; i < worldSize; i++) {
        A_offsets->data()[i] += A_offsets->data()[i-1] + buffWithSizes.data()[i-1].size;
        A_sizes->data()[i] = buffWithSizes.data()[i].size;    
    }
}