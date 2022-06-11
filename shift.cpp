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

// usunieta
// typedef struct shift {
//     int64 i;
//     int64 B_i_h;
// } Shift_data;








void prepareDataForShiftSent(vector<int64>* B, 
                             vector<int64>* None1, 
                             int64 newNodeSize, 
                             int64 nodeSize,
                             int64 dataSize,
                             int64 h,
                             vector<vector<TwoInts64>>* dataForPartitions,
                             int rank,
                             int worldSize) {

    int64 curr_i;
    int64 target_i;
    int targetNode;
    int currNode;

    for (int i = 0; i < B->size(); i++) {
        curr_i = rank * newNodeSize + i;
        target_i = curr_i - h;
        currNode = min((int) (curr_i / newNodeSize), worldSize-1);
        targetNode = min((int) (target_i / newNodeSize), worldSize-1);

        if (target_i >= 0) {
            TwoInts64 data;
            data.i1 = target_i;
            data.i2 = B->data()[i];
            dataForPartitions->data()[targetNode].push_back(data);
        }
        if (curr_i + h > dataSize-1) {
            TwoInts64 data;
            data.i1 = curr_i;
            data.i2 = 0;
            dataForPartitions->data()[currNode].push_back(data);
        }
    }
}


void shift_by_h(vector<int64>** B, 
                vector<int64>** B_new, 
                vector<int64>* SA,
                int64 h,
                int rank, 
                int worldSize) {

    do_sending_operation(B, 
                         B_new, 
                         SA,
                         h, 
                         rank, 
                         worldSize,
                         &prepareDataForShiftSent);
}