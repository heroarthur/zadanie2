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

#include<bits/stdc++.h>

#define root 0


using namespace std;


void prepareDataForShiftSent(vector<int64>* B, 
                             vector<int64>* None1, 
                             int64 newNodeSize, 
                             int64 nodeSize,
                             int64 dataSize,
                             int64 h,
                             vector<vector<TwoInts64>>* dataForPartitions,
                             int rank,
                             int worldSize) {

    // const int64 startIndex = rank * newNodeSize;

    // cout<<"new node size "<<newNodeSize<<endl;

    int64 normalThreadSize;
    int normalVectorUpdateNumber;

    int thread_num;
    int64 threadSize;
    int64 offset;
    int64 threadOffset;
    TwoInts64 data;
    vector<vector<TwoInts64>> localDataForPartitions;
    int curr_i;
    int target_i;
    int currNode;
    int targetNode;
    int data_size_minus_one = dataSize-1;
    int lastNodeIndex = worldSize-1;

    int index;
    int updateVectorsNumber;
    int vector_offset;

    int64 loopEnd;

    assert(B->size() == nodeSize);

    #pragma omp parallel private(loopEnd, normalThreadSize, normalVectorUpdateNumber, curr_i, target_i, currNode, targetNode, thread_num, threadSize, offset, threadOffset, data, localDataForPartitions, index, updateVectorsNumber, vector_offset) num_threads(THREADS_NUM)
    {           

        thread_num = omp_get_thread_num();
        normalThreadSize = ceil(nodeSize / (double) THREADS_NUM);
        normalVectorUpdateNumber = ceil(worldSize / (double) THREADS_NUM);
        threadSize = minInt64(normalThreadSize, maxInt64(0, nodeSize - thread_num * normalThreadSize));

        localDataForPartitions.resize(worldSize);

        offset = rank * newNodeSize;
        threadOffset = thread_num * normalThreadSize;

        // cout<<"zakresy " <<newNodeSize<<" rank "<<rank<<" "<<thread_num<<" start "<<offset<<" koniec "<<offset + threadSize<<" normal thread size "<<threadSize<<endl;
        loopEnd = minInt64(threadOffset + threadSize, nodeSize);

        // cout<<"start koniec "<<offset<<" "<<offset + threadSize<<endl;
        for (int64 i = threadOffset; i < loopEnd; i++) {
            curr_i = offset + i;
            target_i = curr_i - h;

            currNode = min((int) (curr_i / newNodeSize), lastNodeIndex);
            targetNode = min((int) (target_i / newNodeSize), lastNodeIndex);

            if (target_i >= 0) {
                data.i1 = target_i;
                data.i2 = B->data()[i];
                localDataForPartitions.data()[targetNode].push_back(data);
            }
            if (curr_i + h > data_size_minus_one) { //todo zrobic else if
                data.i1 = curr_i;
                data.i2 = 0;
                localDataForPartitions.data()[currNode].push_back(data);
            }
        }

        index = thread_num;

        for (int i = 0; i < THREADS_NUM; i++) {
            index = (index + i) % THREADS_NUM;

            vector_offset = index * normalVectorUpdateNumber;
            updateVectorsNumber = minInt64(normalVectorUpdateNumber, maxInt64(0, worldSize - index * normalVectorUpdateNumber));
            for (int v_index = vector_offset; v_index < vector_offset+updateVectorsNumber; v_index++)
            {
                dataForPartitions->data()[v_index].insert(dataForPartitions->data()[v_index].end(), 
                                                        localDataForPartitions.data()[v_index].begin(), 
                                                        localDataForPartitions.data()[v_index].end());
            }

            #pragma omp barrier
        }
    }
}


void shift_by_h(vector<int64>** B, 
                vector<int64>** B_new, 
                vector<int64>* SA,
                vector<int64>* SA_second_pointer,
                HelpingVectorsSendingOperations* helpVectors,
                int64 h,
                int rank, 
                int worldSize) {

    do_sending_operation(*B, 
                         *B_new, 
                         SA,
                         NULL,
                         helpVectors,
                         h, 
                         false,
                         rank, 
                         worldSize,
                         &prepareDataForShiftSent);

}