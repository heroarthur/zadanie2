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

#ifndef   auxiliary
#define   auxiliary
    #include "auxiliary.cpp"
#endif

#include<bits/stdc++.h>

#define root 0


using namespace std;




void prepareDataForReorderSent(vector<int64>* B, 
                               vector<int64>* SA, 
                               int64 newNodeSize, 
                               int64 nodeSize,
                               int64 dataSize,
                               int64 h,
                               vector<vector<TwoInts64>>* dataForPartitions,
                               int rank,
                               int worldSize) {
    
    #pragma omp parallel for
    for (int thread = 0; thread < worldSize; thread++) {
        for (int i = 0; i < nodeSize; i++) {
            if (thread == getNodeToSend(SA->data()[i], newNodeSize)) {
                TwoInts64 data;
                data.i1 = SA->data()[i];
                data.i2 = B->data()[i];
                dataForPartitions->data()[thread].push_back(data);
            }
        }
    }
}


void reorder_and_rebalance(vector<int64>** B, 
                           vector<int64>** B_new, 
                           vector<int64>* SA,
                           vector<vector<TwoInts64>>* dataForPartitions,
                           HelpingVectors* helpVectors,
                           int rank, 
                           int worldSize) {

    do_sending_operation(*B, 
                         *B_new, 
                         SA,
                         dataForPartitions,
                         helpVectors,
                         EMPTY_HELP_PARAM, 
                         rank, 
                         worldSize,
                         &prepareDataForReorderSent);
    
    switchPointersInt64(B, B_new);
}



