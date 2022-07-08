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


#ifndef   auxiliary
#define   auxiliary
    #include "auxiliary.cpp"
#endif

#include <bits/stdc++.h>

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
    int thread;
    TwoInts64 data;
    // #pragma omp parallel for private(data)
    // for (int thread = 0; thread < worldSize; thread++) {
        #pragma omp parallel for private(thread, data)
        for (int i = 0; i < nodeSize; i++) {
            thread = getNodeToSend(SA->data()[i], newNodeSize);
            // if (thread == getNodeToSend(SA->data()[i], newNodeSize)) {
            data.i1 = SA->data()[i];
            data.i2 = B->data()[i];
            #pragma omp critical 
            {
                dataForPartitions->data()[thread].push_back(data);
            }
            // }
        }
    // }
}


void reorder_and_rebalance(vector<int64>* B, 
                           vector<int64>* B_new, 
                           vector<int64>* SA,
                           vector<int64>* SA_second_pointer,
                           HelpingVectorsSendingOperations* helpVectors,
                           int rank, 
                           int worldSize) {

    do_sending_operation(B, 
                         B_new, 
                         SA,
                         SA_second_pointer,
                         helpVectors,
                         EMPTY_HELP_PARAM,
                         true, 
                         rank, 
                         worldSize,
                         &prepareDataForReorderSent);
    
}



